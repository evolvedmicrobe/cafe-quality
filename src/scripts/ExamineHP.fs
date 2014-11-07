module HPLooker
(* Script to deep dive on the 5 homopolymer error seen so frequently here, and to output covariates related to those with and without errors. *)
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


// The correct and most commonly incorrect template version
//let correct =   "TATACGGGGGGGGGGCACATCA"
//let incorrect = "TATACGGGGGGGGGCACATCA"
let correct =   "CACTGGGGGATAC"
let incorrect = "CACTGGGGATAC"


let correctTemplate = new Sequence(DnaAlphabet.Instance, correct ) //Sequence starts after forty characters, is 22 long
let deletedTemplate = new Sequence(DnaAlphabet.Instance, incorrect)
let correct_rc = correctTemplate.GetReverseComplementedSequence().ConvertToString()
let incorrect_rc = deletedTemplate.GetReverseComplementedSequence().ConvertToString()

// The reference everything is aligned to.
let hpSeq = new Sequence(DnaAlphabet.Instance, "CCCGGGGATCCTCTAGAATGCTCATACACTGGGGGATACATATACGGGGGGGGGGCACATCATCTAGACAGACGACTTTTTTTTTTCGAGCGCAGCTTTTTGAGCGACGCACAAGCTTGCTGAGGACTAGTAGCTTC")
let hpRef = new  Reference(hpSeq);

// Calculate the size of the indel in the 5 bp homopolymer
let DeletionSize (sub:Sequence) =
   let alns = hpRef.AlignSequence (sub) |> Seq.toArray
   if alns.Length  = 0 then "NaN" else
       let best = alns.[0]       
       let variants = VariantCaller.CallVariants (best, hpRef.RefSeq)
       let bad = variants |> Seq.where (fun u -> u.StartPosition = 29) |> Seq.toArray
       if bad.Length = 0 then "0" else
           if bad.[0].Type = VariantType.SNP then "SNP" else
               let mut = bad.[0] :?> IndelVariant
               match mut.InsertionOrDeletion with
                   | IndelType.Deletion -> "-" + mut.InsertedOrDeletedBases.Length.ToString()
                   | IndelType.Insertion -> mut.InsertedOrDeletedBases.Length.ToString()
                   | _ -> failwith "bad"

let cntMergeSpikes (read: ReadFromZMW)  =
    let bp = match read.ReverseComplementedOriginally with
             | true -> 'G'
             | false -> 'C'
    // Now to find where the homopolymer is inside sequence
    (read.BaseCalls, read.MergeQV) ||> Seq.map2 (fun x y -> x = bp && y > (byte 70)) |> Seq.where id |> Seq.length

(* Take a subread and only extract out the bit that matches the 5 bp homopolymer plus a window around the sides *)
let getHPSection (parentRead : ReadFromZMW) (subRead : CCSSubRead) = 
   let toAlign = new Sequence(DnaAlphabet.Instance, subRead.Seq,false)
   let mutable alns = hpRef.AlignSequence(toAlign)
   if alns.Count = 0 then None else
       let mutable top = alns.[0]
       //Flip the read if necessary
       let revComp = top.SecondSequence.Metadata.ContainsKey("+isReversed")
       let mutable zmwSection = parentRead.GetSubRead(subRead)
       //if revComp then 
        //   zmwSection <- zmwSection.GetReverseComplement()
        //   top <- hpRef.AlignSequence( (toAlign.GetReverseComplementedSequence() :?> Sequence )).[0];
       //See if it overlaps with 5 bp homopolymer and output if so.
       let start_pos = top.FindQueryPositionCorrespondingtoReferencePosition(26) 
       let end_pos = top.FindQueryPositionCorrespondingtoReferencePosition(26+11) 
       if start_pos.HasValue && end_pos.HasValue then
           let start = int start_pos.Value
           let endi = int end_pos.Value
           let len = endi - start + 1
           let subSection = match revComp with
                            | false -> zmwSection.GetSubSection(start, len)
                            | true -> zmwSection.GetSubSection((int subRead.Seq.Length) - len - start, len)
           subSection.ReverseComplementedOriginally <- revComp
           subSection.SpikeMergeQVCount <- cntMergeSpikes subSection
           subSection.OriginalSubReadLength <- subRead.Seq.Length
           subSection.HomopolymerLength <- (DeletionSize (new Sequence(DnaAlphabet.Instance, subRead.Seq, false)))
           subSection.RQ <- subRead.RQ
           subSection.Zmw <- subRead.ParentZMW
           Some(subSection) else
           None

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

let outFile = csv_writer("/Users/nigel/git/cafe-quality/data/homopolymerDeepDive5bpLong.csv")




//Code to calculate various likelihoods.
let qConfig = TemplateRegion.GetQuiverConfig()
let viterbi = new SparseSseQvReadScorer(qConfig)
let sumproduct = new SparseSseQvSumProductReadScorer(qConfig)
let makeQuiverRead (read: ReadFromZMW) =
   let qs = new QvSequenceFeatures(read.BaseCalls,read.InsertionQV,read.SubstitutionQV, read.DeletionQV, read.DeletionTag, read.MergeQV)
   new Read(qs, "tmp", "P6-C4")


//Map all the reads by taking all of their subreads.  
let ProcessRead (read:CCSRead) = 
   let mutable badRead = new ReadFromZMW()
   let del_size = DeletionSize read.Seq
   if del_size <> "NaN" then
       read.ZMW <- CCSRead.ParentExperiment.GetZMWforRead(read)
       let readzmw = new ReadFromZMW(read.ZMW.UnsafeFullRead)
       let procesSubRead = getHPSection readzmw
       let hps = read.SubReads |> Seq.truncate 126 |> Seq.map procesSubRead |> Seq.toArray
       let mutable subRead = 0
       for region in hps do
           subRead <- subRead + 1
           if region.IsSome then
               let v = region.Value
               v.ConsensusIndelSize <- del_size
               let qread = makeQuiverRead v
               v.SubReadNumber <- subRead
               v.HPSectionLength <-  v.BaseCalls.Length
               try
                   if not v.ReverseComplementedOriginally then
                       v.NoErrorViterbiScore <- viterbi.Score(correct, qread)
                       v.OneDeletionErrorViterbiScore <- viterbi.Score(incorrect, qread)
                       v.NoErrorSumProductScore <- sumproduct.Score(correct, qread)
                       v.OneDeletionSumProductScore <- sumproduct.Score(incorrect,qread)
                    else
                       v.NoErrorViterbiScore <- viterbi.Score(correct_rc, qread)
                       v.OneDeletionErrorViterbiScore <- viterbi.Score(incorrect_rc, qread)
                       v.NoErrorSumProductScore <- sumproduct.Score(correct_rc, qread)
                       v.OneDeletionSumProductScore <- sumproduct.Score(incorrect_rc,qread)
                   
                with
                    | _ -> ()
               let outData = OutputHelper.CalculateDataLines(v)
               outFile.write outData
           else
               badRead.Zmw <- read.ZMWnumber
               badRead.SubReadNumber <- subRead
               badRead.ConsensusIndelSize <- del_size
               let outData = OutputHelper.CalculateDataLines(badRead)
               outFile.write outData


[<EntryPoint>]
let main args =
   let sw = System.Diagnostics.Stopwatch.StartNew()
   LoadZMWs.ccs_data.CCSReads  |> Seq.where (fun z -> z.AssignedReference <> null && z.AssignedReference.RefSeq.ID = "HP.V1.02" && z.SubReads.Count > 126) 
                               |> Seq.truncate 5000
                               |> Seq.toArray
                               //|> Array.iter (fun j -> ProcessRead j)
                               |> Array.Parallel.iter (fun j -> ProcessRead j)
   sw.Stop()
   outFile.Close
   printfn "%f" sw.Elapsed.TotalSeconds
   // Now to ouptu
   Console.WriteLine("Success");
   0
