namespace deletem
//(* A quick script to grab several reads from the HP experiment and periodically output a subsample of QV values *)
//module QVFinder
//open System
//open System.IO
//open VariantCaller
//open Bio
//open Bio.Algorithms.Alignment
//open PacBio.Data
//open PacBio.Utils
//open ConsensusCore
//open LoadZMWs
//open System.Linq
//open PacBio.Consensus
//
//// Subsample down to this amount
//let percentageReadsToUse = 0.05
//let rand = new Random()
//
//let getSection (parentRead : ReadFromZMW) (subRead : CCSSubRead) = 
//    if rand.NextDouble() < percentageReadsToUse then None else
//        let zmwSection = parentRead.GetSubRead(subRead)
//        zmwSection.Zmw = parentRead.Zmw
//        Some(zmwSection)
//
//// Output Relevant Data
//type csv_writer (fname:string) =
//    let sw = new StreamWriter(fname)
//    let header = "Zmw,Base,MergeQV,DeletionQV,DeletionTag,SubsQv,InsQV"
//    do  sw.WriteLine(header)
//    member this.write data =
//            let toW = data
//                        |> Seq.map (fun x -> x.ToString())
//                        |> String.concat ","
//            lock sw ( fun () -> sw.WriteLine (toW))
//    member this.writeReadAtBase (read : ReadFromZMW) (i: int) = 
//                let data : obj[] = [| read.Zmw; read.BaseCalls.[i]; read.MergeQV.[i]; read.DeletionQV.[i]; read.DeletionTag.[i];
//                                      read.SubstitutionQV.[i];read.InsertionQV.[i]; |]
//                this.write data
//                
//    member this.Close = sw.Close()
//
//let outFile = csv_writer("/Users/nigel/git/cafe-quality/data/QVvalues.csv")
//
//
//// Map all the reads by taking all of their subreads.  
//let ProcessRead (read:CCSRead) = 
//    read.ZMW <- CCSRead.ParentExperiment.GetZMWforRead(read)
//    let readzmw = new ReadFromZMW(read.ZMW.UnsafeFullRead)
//    let procesSubRead = getSection readzmw
//    let subs = read.SubReads |> Seq.map procesSubRead 
//                            |> Seq.choose id 
//    let processSub (sub : ReadFromZMW) = 
//        for i in 0 .. (sub.BaseCalls.Length - 1) do
//            outFile.writeReadAtBase sub i
//    subs |> Seq.iter processSub
//
//[<EntryPoint>]
//let main args =
//    let sw = System.Diagnostics.Stopwatch.StartNew()
//    LoadZMWs.ccs_data.CCSReads  |> Seq.where (fun z -> z.AssignedReference <> null && z.AssignedReference.RefSeq.ID = "HP.V1.02") 
//                                |> Seq.truncate 2000
//                                |> Seq.toArray
//                                |> Array.Parallel.iter (fun j -> ProcessRead j)
//    sw.Stop()
//    outFile.Close
//    printfn "%f" sw.Elapsed.TotalSeconds
//    // Now to ouptu
//    Console.WriteLine("Success");
//    0
