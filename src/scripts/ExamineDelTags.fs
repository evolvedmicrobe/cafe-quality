module DelTagInvestigation
(* Script to examine how frequently DelTags are used and how often they are correct. *)
open System
open System.IO
open VariantCaller
open Bio
open Bio.Algorithms.Alignment
open PacBio.Data
open PacBio.Utils
open ConsensusCore
open LoadZMWs
open System.Linq
open PacBio.Consensus
open System.Collections.Generic
open Microsoft.FSharp.Collections

type DelTagReport = { Length: int; DelTags: int; CorrectDelTags: int}



let countCorrectDelTags (aln : PairwiseAlignedSequence) (subSection : ReadFromZMW) =
    // Presumably at this stage the alignment and the ReadFromZMW match.
    // Verify the alignment and get an array to translate between them.
    let mutable seq_pos = -1
    let translationArray = Array.init ((int) aln.SecondSequence.Count) (fun i -> -999)
    for aln_pos in 0 .. ((int)aln.SecondSequence.Count - 1) do
        let bp = aln.SecondSequence.[(int64)aln_pos]
        if bp <> (byte)'-' then 
            seq_pos <- seq_pos + 1
            translationArray.[aln_pos] <- seq_pos
            if bp <> (byte)subSection.BaseCalls.[seq_pos] then 
                failwith "Alignment didn't match"
        
    let tag_positions = seq { for i in 0 .. (subSection.DeletionTag.Length - 1) do if subSection.DeletionTag.[i] <> (byte)'N' then yield i} |> Seq.toList
    // Now for each DelTag I need to determine if it is valid.
    let mutable correctTags = 0
    for seq_position in tag_positions do
        if seq_position <> 0 then
            // Check for the simple case of a deletion right before this base.
            let del_tag = subSection.DeletionTag.[seq_position]
            let aln_position = (Array.findIndex (fun z -> z = seq_position) translationArray) - 1
            if aln.SecondSequence.[(int64)aln_position] = (byte)'-' && del_tag = aln.FirstSequence.[(int64)aln_position] then 
                correctTags <- correctTags + 1 
            else
                // Failing that, we might have a more complicated case possible due to the alignment being skewed, this generally sucks to handle
                // As a simple rule, I will only try to handle simple homopolymer cases, if I can shift the base over, great
                // This does not strictly account for gap open penalties, and all cases where the deletion tag can't be proven to be wrong
                let c_pos = aln_position
                let first_gap = Seq.tryFindIndex (fun j -> j = (byte) 'N') (aln.SecondSequence.Take( (aln_position + 1)))
                if first_gap.IsSome then 
                    let gap_pos = first_gap.Value
                    // Only fixing single gaps 
                    if ((gap_pos - 1) > 0 ) && aln.SecondSequence.[((int64)gap_pos - 1L)] <> (byte)'N' then 
                        // Single basepair gap, can I move it down to one basepair before the deletion without a change in score?
                        // criteria has to be all bases in reference are the same
                        let bases_present = seq { for i in gap_pos .. aln_position do yield aln.FirstSequence.[(int64)i] } |> Seq.distinct |> Seq.length
                        if bases_present = 0 then
                            // check if it's a match, assume we can move over below
                            let deleted = aln.FirstSequence.[(int64)aln_position]
                            if deleted = del_tag then correctTags <- correctTags + 1
    {Length = (int)subSection.BaseCalls.Length; DelTags = tag_positions.Length; CorrectDelTags = correctTags}
    
    


(* Take a subread and only extract out the bit that matches the 5 bp homopolymer plus a window around the sides *)
let getTagReport (parentRead : ReadFromZMW) (subRead : CCSSubRead) =
   Console.WriteLine( subRead.ParentZMW)
   let toAlign = new Sequence(DnaAlphabet.Instance, subRead.Seq, false)
   let mutable alns = parentRead.AssignedReference.AlignSequence(toAlign)
   if alns.Count = 0 then 
       None else
       //Console.WriteLine(alns.[0].ToString())
       let mutable top = alns.[0]
       //Flip the read if necessary
       let revComp = top.SecondSequence.Metadata.ContainsKey("+isReversed")
       let mutable zmwSection = parentRead.GetSubRead(subRead)
       if revComp then 
          zmwSection <- zmwSection.GetReverseComplement()
          top <- parentRead.AssignedReference.AlignSequence( (toAlign.GetReverseComplementedSequence() :?> Sequence )).[0];
        
       // Get the subsection corresponding to the alignment
       let start = int top.SecondSequenceStart.Value
       let len = top.SecondSequence |> Seq.where (fun u -> u <> (byte)'-') |> Seq.length
       let subSection = zmwSection.GetSubSection(start, len)

       let del_rpt = countCorrectDelTags top subSection
       subSection.AssignedReference <- parentRead.AssignedReference
       subSection.ReverseComplementedOriginally <- revComp
       subSection.OriginalSubReadLength <- subRead.Seq.Length
       subSection.RQ <- subRead.RQ
       subSection.Zmw <- subRead.ParentZMW
       subSection.CountCorrectDelTags <- del_rpt.CorrectDelTags
       subSection.CountDelTags <- del_rpt.DelTags
       subSection.AlignedLength <- del_rpt.Length
       Some(subSection)

//Output Relevant Data
type csv_writer (fname:string) =
   let sw = new StreamWriter(fname)
   let header = ReadFromZMW.CSVHeaderFields() |> String.concat ","
   do  sw.WriteLine(header)
   member this.write data =
           let toW = data
                       |> Seq.map (fun x -> x.ToString())
                       |> String.concat ","
           lock sw ( fun () -> sw.WriteLine (toW))
   
   member this.Close = sw.Close()

let outFile = csv_writer("/Users/nigel/git/cafe-quality/data/DelTagInvestigation.csv")
let mutable totalFailures = 0
//Map all the reads by taking all of their subreads.  
let ProcessRead (read:CCSRead)  = 
   try
       let mutable badRead = new ReadFromZMW()
       read.ZMW <- CCSRead.ParentExperiment.GetZMWforRead(read)
       let readzmw = new ReadFromZMW(read.ZMW.UnsafeFullRead)
       readzmw.AssignedReference <- read.AssignedReference
       let procesSubRead = getTagReport readzmw
       let hps = read.SubReads |> Seq.map procesSubRead |> Seq.toArray
       let mutable subRead = 0
       for region in hps do
           subRead <- subRead + 1
           if region.IsSome then
               let v = region.Value
               let outData = OutputHelper.CalculateDataLines(v)
               outFile.write outData
           else
               badRead.Zmw <- read.ZMWnumber
               badRead.SubReadNumber <- subRead
               let outData = OutputHelper.CalculateDataLines(badRead)
               outFile.write outData
    with
        | _ -> totalFailures <- totalFailures + 1


let main =
   let sw = System.Diagnostics.Stopwatch.StartNew()
   let refs = LoadZMWs.ccs_data.References
   LoadZMWs.ccs_data.CCSReads  |> Seq.where (fun z -> z.AssignedReference <> null && z.SubReads.Count > 1 && z.AssignedReference.RefSeq.ID <> "SmrtBellSequence") 
                               |> Seq.truncate 1000
                               |> Seq.toArray
                               |> Array.iter (fun j -> ProcessRead j)
                               //|> Array.Parallel.iter (fun j -> ProcessRead j)
   sw.Stop()
   outFile.Close
   printfn "%f" sw.Elapsed.TotalSeconds
   // Now to ouptu
   Console.WriteLine("Success");
   Console.WriteLine("Failures: "+totalFailures.ToString())
   0
