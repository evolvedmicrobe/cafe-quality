module DumpForLisa
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
let all4Mer = new Sequence(DnaAlphabet.Instance,"GAAGCTACTAGTCCTCAGCAAGCTTGGATCAACACTAAATTATGCGAGCGCGGATAGATTCGCTCTGCAATCCGTCCCACGATGTGTTGGTATATTGACTGGGAGAGTTCATAACCCTATCTCCTTGCTACGCCTCGAACGGGGTTTAAGGCCCCGGCTTCCATTTTGTCGGTCTAGTAATGGCGTGAAGACCGACGTTAGGGCAGAATACATGAGGTGGACAGGAAACTCAGCCAAAAGTGCCGCACCAGTCACAAGCTGTAGCATCGTACTTTCTTACCTGATGCATTCTAGAGGATCCCCGGG")
let all4MerRef = new  Reference(all4Mer);

(* Take a subread and extract the template and portion of the read aligned to it *)
let getOutputRead (parentRead : ReadFromZMW) (subRead : CCSSubRead) = 
   try
       let toAlign = new Sequence(DnaAlphabet.Instance, subRead.Seq,false)
       let mutable alns = all4MerRef.AlignSequence(toAlign)
       if alns.Count <> 1 then 
           None else
           let mutable top = alns.[0]
           let tmp_tpl = top.FirstSequence.Where( (fun z -> z <> (byte '-'))) |> Seq.map (fun u -> char u) |> Seq.toArray
           let tmp_rd = top.SecondSequence.Where( (fun z -> z <> (byte '-'))) |> Seq.map (fun u -> char u) |> Seq.toArray
           let mutable tpl = new String(tmp_tpl)
           let mutable rd  = new String(tmp_rd)
           //Flip the read if necessary
           let revComp = top.SecondSequence.Metadata.ContainsKey("+isReversed")
           if revComp then 
                tpl <- (new Sequence(NoGapDnaAlphabet.Instance, tpl)).GetReverseComplementedSequence().ConvertToString()
                rd <- (new Sequence(NoGapDnaAlphabet.Instance, rd)).GetReverseComplementedSequence().ConvertToString()
           Some((tpl, rd))
    with
        | _ -> None
        
     

//Output Relevant Data
type csv_writer (fname:string) =
   let sw = new StreamWriter(fname)
   let header = "ZMW,Subread,Template,Read"
   do  sw.WriteLine(header)

   member this.write (zmw:Object) (i:Object) (tpl:Object) (read:Object) =
           let toW = [|zmw; i; tpl; read|]
                       |> Array.map (fun x -> x.ToString())
                       |> String.concat ","
           lock sw ( fun () -> sw.WriteLine (toW))
   
   member this.Close = sw.Close()

let outFile = csv_writer("/Users/nigel/git/cafe-quality/data/TemplateReadPairs.csv")


//Map all the reads by taking all of their subreads.  
let ProcessRead (read:CCSRead) = 
   let mutable badRead = new ReadFromZMW()
   read.ZMW <- CCSRead.ParentExperiment.GetZMWforRead(read)
   let readzmw = new ReadFromZMW(read.ZMW.UnsafeFullRead)
   let procesSubRead = getOutputRead readzmw
   let subs = read.SubReads |> Seq.truncate 126 |> Seq.map procesSubRead |> Seq.toArray
   let mutable subRead = 0
   for rd_tpl in subs do
       subRead <- subRead + 1
       if rd_tpl.IsSome then outFile.write read.ZMWnumber subRead (fst rd_tpl.Value) (snd rd_tpl.Value)
   1


let main =
   let sw = System.Diagnostics.Stopwatch.StartNew()
   let reads = LoadZMWs.ccs_data.CCSReads  |> Seq.where (fun z -> z.AssignedReference <> null && z.AssignedReference.RefSeq.ID = "ALL4MER.V2.01") |> Seq.toList
   Console.WriteLine("Count " + reads.Length.ToString())
   let tot = reads |> Seq.map (fun j -> ProcessRead j) |> Seq.sum
   Console.WriteLine(tot)
   sw.Stop()
   outFile.Close
   printfn "%f" sw.Elapsed.TotalSeconds
   // Now to ouptu
   Console.WriteLine("Success");
   0
