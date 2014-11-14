//namespace deleteme
//
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



// The reference everything is aligned to.
let hpSeq = new Sequence(DnaAlphabet.Instance, "CCCGGGGATCCTCTAGAATGCTCATACACTGGGGGATACATATACGGGGGGGGGGCACATCATCTAGACAGACGACTTTTTTTTTTCGAGCGCAGCTTTTTGAGCGACGCACAAGCTTGCTGAGGACTAGTAGCTTC")
let hpRef = new  Reference(hpSeq);

(* Take a subread and only extract out the bit that matches the 5 bp homopolymer plus a window around the sides *)
let getHPSection (parentRead : ReadFromZMW) (subRead : CCSSubRead) = 
    let toAlign = new Sequence(DnaAlphabet.Instance, subRead.Seq,false)
    let mutable alns = hpRef.AlignSequence(toAlign)
    if alns.Count = 0 then None else
        let mutable top = alns.[0]
        // Flip the read if necessary
        let revComp = top.SecondSequence.Metadata.ContainsKey("+isReversed")
        if revComp then None else
            let mutable zmwSection = parentRead.GetSubRead(subRead)
            // See if it overlaps with 5 bp homopolymer and output if so.
            let start_pos = top.FindQueryPositionCorrespondingtoReferencePosition(40)
            let end_pos = top.FindQueryPositionCorrespondingtoReferencePosition(40+22)
            if start_pos.HasValue && end_pos.HasValue then
                let start = (int start_pos.Value) + subRead.Start
                let endi = (int end_pos.Value) + subRead.Start
                Some((start, endi)) else
                None

// Output Relevant Data
type csv_writer (fname:string) =
    let sw = new StreamWriter(fname)
    let header = "Movie,HoleNumber,StartBase,StopBase"
    do  sw.WriteLine(header)
    member this.write data =
            let toW = data
                        |> Seq.map (fun x -> x.ToString())
                        |> String.concat ","
            lock sw ( fun () -> sw.WriteLine (toW))    
    member this.Close = sw.Close()

let outFile = csv_writer("/Users/nigel/git/cafe-quality/data/HP_Sections.csv")

// Map all the reads by taking all of their subreads.  
let ProcessRead (read:CCSRead) = 
    let mutable badRead = new ReadFromZMW()
    read.ZMW <- CCSRead.ParentExperiment.GetZMWforRead(read)
    let readzmw = new ReadFromZMW(read.ZMW.UnsafeFullRead)
    let procesSubRead = getHPSection readzmw
    let hps = read.SubReads |> Seq.truncate 126 |> Seq.map procesSubRead |> Seq.toArray
    let mutable subRead = 0
    for region in hps do
        subRead <- subRead + 1
        if region.IsSome then
            let v = region.Value
            let outData = [| read.Movie; read.ZMWnumber.ToString(); (fst v).ToString(); (snd v).ToString() |]
            outFile.write outData



[<EntryPoint>]
let main args =
    let sw = System.Diagnostics.Stopwatch.StartNew()
    LoadZMWs.ccs_data.CCSReads  |> Seq.where (fun z -> z.AssignedReference <> null && z.AssignedReference.RefSeq.ID = "HP.V1.02" && z.SubReads.Count > 126) 
                                |> Seq.truncate 200
                                |> Seq.toArray
                                |> Array.Parallel.iter (fun j -> ProcessRead j)
    sw.Stop()
    outFile.Close
    printfn "%f" sw.Elapsed.TotalSeconds
    // Now to ouptu
    Console.WriteLine("Success");
    0
