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
let ref_str = "CCCGGGGATCCTCTAGAATGCTCATACACTGGGGGATACATATACGGGGGGGGGGCACATCATCTAGACAGACGACTTTTTTTTTTCGAGCGCAGCTTTTTGAGCGACGCACAAGCTTGCTGAGGACTAGTAGCTTC"
let correct =   "CACTGGGGGATAC"
let incorrect = "CACTGGGGATAC"
let hp_char = 'G'
let start_index = 29
//let correct =   "CAGCTTTTTGAGC"
//let incorrect = "CAGCTTTTGAGC"

let correctTemplate = new Sequence(DnaAlphabet.Instance, correct ) //Sequence starts after forty characters, is 22 long
let deletedTemplate = new Sequence(DnaAlphabet.Instance, incorrect)
let correct_rc = correctTemplate.GetReverseComplementedSequence().ConvertToString()
let incorrect_rc = deletedTemplate.GetReverseComplementedSequence().ConvertToString()

// The reference everything is aligned to.
let hpSeq = new Sequence(DnaAlphabet.Instance,         "CCCGGGGATCCTCTAGAATGCTCATACACTGGGGGATACATATACGGGGGGGGGGCACATCATCTAGACAGACGACTTTTTTTTTTCGAGCGCAGCTTTTTGAGCGACGCACAAGCTTGCTGAGGACTAGTAGCTTC")
let deletionHPSeq = new Sequence(DnaAlphabet.Instance, "CCCGGGGATCCTCTAGAATGCTCATACACTGGGGATACATATACGGGGGGGGGGCACATCATCTAGACAGACGACTTTTTTTTTTCGAGCGCAGCTTTTTGAGCGACGCACAAGCTTGCTGAGGACTAGTAGCTTC")

let fullCorrect = hpSeq.ConvertToString()
let fullCorrect_rc = hpSeq.GetReverseComplementedSequence().ConvertToString()
let full_deletion = deletionHPSeq.ConvertToString()
let full_deletion_rc = deletionHPSeq.GetReverseComplementedSequence().ConvertToString()

let hpRef = new  Reference(hpSeq);

let convertToProb x = Math.Pow(-x / 10.0,10.0)

(* Calculate the size of the indel in the 5 bp homopolymer
   I am doing this two ways, one by aligning and calling variants, 
   and a second by counting the homopolymers in the region, to account for SNPs within
   which is much mor common than I would have thought.
*)
let DeletionSize (sub:Sequence) =
   let alns = hpRef.AlignSequence (sub) |> Seq.toArray
   if alns.Length  = 0 then "NaN" else
       let best = alns.[0]       
       // First get HP length determined by alignment
       let variants = VariantCaller.CallVariants (best, hpRef.RefSeq)
       let bad = variants |> Seq.where (fun u -> u.StartPosition = start_index ) |> Seq.toArray
       if bad.Length = 0 then "0" else
           if bad.[0].Type = VariantType.SNP then "SNP" else
               let mut = bad.[0] :?> IndelVariant
               match mut.InsertionOrDeletion with
                   | IndelType.Deletion -> "-" + mut.InsertedOrDeletedBases.Length.ToString()
                   | IndelType.Insertion -> mut.InsertedOrDeletedBases.Length.ToString()
                   | _ -> failwith "bad"


(* To avoid alignment issues, I am going to just grab the region, look for the longest continuous stretch of 
    the homopolymer base, and count the number of merge spikes, deletion tags and its length *)
type homoPolymerReport = { Length: int; MergeSpikes: int; DeletionTags: int; ErrorProb: double}

let cntHPEvents (read: ReadFromZMW)  =
    // First task, where is the longest homopolymer stretch ?
    let bp = match read.ReverseComplementedOriginally with
             | true -> 'C'
             | false -> 'G'
    let i = ref 0
    // Count the start and length of each position for the homopolymer base
    let groups = seq { while !i < read.BaseCalls.Length do
                        if read.BaseCalls.[!i] = bp then
                            let start = !i
                            let cnt = ref 1
                            i := !i + 1
                            while !i < read.BaseCalls.Length && read.BaseCalls.[!i] = bp do
                                cnt := !cnt + 1
                                i := !i + 1
                            yield (start, !cnt) else i := !i + 1 } |> Seq.toList
    if groups.Length = 0 then
        {Length = -99; MergeSpikes = -999; DeletionTags = -999} else
        let (startHP, Length) = groups |> Seq.maxBy snd
        // Now to find where the homopolymer is inside sequence
        let hpRange = seq {startHP .. (startHP + Length - 1) } 
        if (Seq.max hpRange) > read.BaseCalls.Length then failwith "Mother fucker"
        let spikeCount =    hpRange |> 
                                Seq.map (fun i -> read.BaseCalls.[i] = bp && read.MergeQV.[i] > (byte 70)) |> 
                                Seq.where id |> Seq.length
        let delTagCount = hpRange |> Seq.where (fun j -> (j+1) < read.BaseCalls.Length && read.DeletionTag.[j+1] = (byte)bp) |> Seq.length
        let expectedErrors = hpRange |> Seq.map (fun x -> convertToProb(read.MergeQV[x]) + convertToProb(read.DeletionQV[x]) 
        {Length = Length; MergeSpikes = spikeCount; DeletionTags = delTagCount}


(* Take a subread and only extract out the bit that matches the 5 bp homopolymer plus a window around the sides *)
let getHPSection (parentRead : ReadFromZMW) (subRead : CCSSubRead) = 
   let toAlign = new Sequence(DnaAlphabet.Instance, subRead.Seq,false)
   let mutable alns = hpRef.AlignSequence(toAlign)
   if alns.Count = 0 then 
       None else
       //Console.Write("Full")
       //Console.WriteLine(alns.[0].ToString())
       let mutable top = alns.[0]
       //Flip the read if necessary
       let revComp = top.SecondSequence.Metadata.ContainsKey("+isReversed")
       let mutable zmwSection = parentRead.GetSubRead(subRead)
       //if revComp then 
        //   zmwSection <- zmwSection.GetReverseComplement()
        //   top <- hpRef.AlignSequence( (toAlign.GetReverseComplementedSequence() :?> Sequence )).[0];
       //See if it overlaps with 5 bp homopolymer and output if so.
       let start_pos = top.FindQueryPositionCorrespondingtoReferencePosition(26) 
       let end_pos = top.FindQueryPositionCorrespondingtoReferencePosition(26 + 11) 
       if start_pos.HasValue && end_pos.HasValue then
           let start = int start_pos.Value
           let endi = int end_pos.Value
           let len = endi - start + 1
           let subSection = match revComp with
                            | false -> zmwSection.GetSubSection(start, len)
                            | true -> zmwSection.GetSubSection((int subRead.Seq.Length) - len - start, len)
           subSection.ReverseComplementedOriginally <- revComp
           let hpEvents = cntHPEvents subSection
           subSection.SpikeMergeQVCount <- hpEvents.MergeSpikes
           subSection.HPDeletionTagNotNCount <- hpEvents.DeletionTags
           subSection.HomopolymerLengthFromCounting <- hpEvents.Length
           subSection.OriginalSubReadLength <- subRead.Seq.Length
           subSection.HomopolymerLengthFromAlignment <- (DeletionSize (new Sequence(DnaAlphabet.Instance, subRead.Seq, false)))
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

let outFile = csv_writer("/Users/nigel/git/cafe-quality/data/homopolymerDeepDiveDiagnostics5bpP6.csv")




//Code to calculate various likelihoods.
let qConfigP6 = TemplateRegion.GetP6C4QuiverConfig()
let viterbiP6 = new SparseSseQvReadScorer(qConfigP6)
let sumproductP6 = new SparseSseQvSumProductReadScorer(qConfigP6)

let qConfigC2 = TemplateRegion.GetC2QuiverConfig()
let viterbiC2 = new SparseSseQvReadScorer(qConfigC2)
let sumproductC2 = new SparseSseQvSumProductReadScorer(qConfigC2)

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
           let cSubRead = read.SubReads.[subRead]
           subRead <- subRead + 1
           let fread = makeQuiverRead (readzmw.GetSubRead(cSubRead))
           if region.IsSome then
               let v = region.Value
               v.ConsensusIndelSize <- del_size
               let qread = makeQuiverRead v
               v.SubReadNumber <- subRead
               v.HPSectionLength <-  v.BaseCalls.Length
               try
                   if not v.ReverseComplementedOriginally then
                       v.NoErrorViterbiScoreP6 <- viterbiP6.Score(correct, qread)
                       v.OneDeletionErrorViterbiScoreP6 <- viterbiP6.Score(incorrect, qread)
                       v.NoErrorSumProductScoreP6 <- sumproductP6.Score(correct, qread)
                       v.OneDeletionSumProductScoreP6 <- sumproductP6.Score(incorrect,qread)
                       v.FullLengthCorrectSumProductScore <- sumproductP6.Score(fullCorrect, fread)
                       v.FullLengthIncorrectSumProductScore <- sumproductP6.Score(full_deletion, fread)
                       v.SummedPulseWidthForHP <- v.BaseCalls |> Seq.mapi (fun i c -> if c = 'G' then (float32)v.IpdInFrames.[i] else 0.0f) |> Seq.sum
                    

                       v.NoErrorViterbiScoreC2 <- viterbiC2.Score(correct, qread)
                       v.OneDeletionErrorViterbiScoreC2 <- viterbiC2.Score(incorrect, qread)
                       v.NoErrorSumProductScoreC2 <- sumproductC2.Score(correct, qread)
                       v.OneDeletionSumProductScoreC2 <- sumproductC2.Score(incorrect,qread)
                  
                    else
                       v.NoErrorViterbiScoreP6 <- viterbiP6.Score(correct_rc, qread)
                       v.OneDeletionErrorViterbiScoreP6 <- viterbiP6.Score(incorrect_rc, qread)
                       v.NoErrorSumProductScoreP6 <- sumproductP6.Score(correct_rc, qread)
                       v.OneDeletionSumProductScoreP6 <- sumproductP6.Score(incorrect_rc,qread)
                       v.FullLengthCorrectSumProductScore <- sumproductP6.Score(fullCorrect_rc, fread)
                       v.FullLengthIncorrectSumProductScore <- sumproductP6.Score(full_deletion_rc, fread)
                       v.SummedPulseWidthForHP <- v.BaseCalls |> Seq.mapi (fun i c -> if c = 'C' then (float32)v.IpdInFrames.[i] else 0.0f) |> Seq.sum
                
                       v.NoErrorViterbiScoreC2 <- viterbiC2.Score(correct_rc, qread)
                       v.OneDeletionErrorViterbiScoreC2 <- viterbiC2.Score(incorrect_rc, qread)
                       v.NoErrorSumProductScoreC2 <- sumproductC2.Score(correct_rc, qread)
                       v.OneDeletionSumProductScoreC2 <- sumproductC2.Score(incorrect_rc,qread)
                     

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


//[<EntryPoint>]
let main args =
   let sw = System.Diagnostics.Stopwatch.StartNew()
   LoadZMWs.ccs_data.CCSReads  |> Seq.where (fun z -> z.AssignedReference <> null && z.AssignedReference.RefSeq.ID = "HP.V1.02" && z.SubReads.Count > 126) 
                               |> Seq.truncate 1000
                               |> Seq.toArray
                               //|> Array.iter (fun j -> ProcessRead j)
                               |> Array.Parallel.iter (fun j -> ProcessRead j)
   sw.Stop()
   outFile.Close
   printfn "%f" sw.Elapsed.TotalSeconds
   // Now to ouptu
   Console.WriteLine("Success");
   0
