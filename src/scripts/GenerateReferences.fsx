 (* 
   This is a file that I used to get the reference sequence data from PacBio.Data.   
   I took the top three peaks from the CCS read length distributions and mapped them back to the reference.  
   There were several minor differences, in particular with the All4mers.   
 *)
#I  "/Users/nigel/git/cafe-quality/lib/" 
#r "PacBio.Utils.dll"
#r "PacBio.HDF.dll"
#r "PacBio.IO.dll"
#r "VariantCaller.dll"
#r "Bio.dll"

open System;
open System.IO
open VariantCaller
open PacBio.HDF
open PacBio.Utils
open PacBio.IO
open Bio.Util
open Bio

let data = 
    let runninMono = Bio.CrossPlatform.Environment.RunningInMono
    let d_direct = if runninMono  then 
                        "/Users/nigel/git/cafe-quality/data/" else
                        @"C:\git\cafe-quality\data\"

    let direc = if runninMono then
                    "/Users/nigel/CCS_P6_C4/Analysis_Results/"
                else
                    @"C:\CCS_P6_C4\Analysis_Results\\"

    let ccsFiles = (new DirectoryInfo(d_direct)).GetFiles() |> 
                    Seq.where (fun h -> h.Name.EndsWith(".ccs.fasta.gz")) |> 
                    Seq.map (fun u-> u.FullName) |>
                    Seq.toList

    let reference = Path.Combine(d_direct, "References.fna")

    let subReads = (new DirectoryInfo(d_direct)).GetFiles() |> 
                    Seq.where (fun h -> h.Name.EndsWith(".subreads.fasta.gz")) |> 
                    Seq.map (fun u-> u.FullName) |>
                    Seq.toList

    let qc_exp = new QualityExperiment.QualityExperiment(null, ccsFiles, subReads, reference)
    Console.WriteLine qc_exp.CCSReads.Count
    qc_exp



let alnLength = data.References |> Seq.map (fun ref -> ref.ReferenceSequence.Count) |> Seq.toArray

let sw = new StreamWriter("/Users/nigel/git/cafe-quality/lib/CCS_Lengths.csv")
sw.WriteLine("ZMW,ReadType,Length,RQ")
let outRead (read:CCSRead) =
   let zmw = read.ZMW.ToString()
   sw.WriteLine(zmw+",0,"+read.Seq.Count.ToString()+",NaN")
   read.SubReads |> Seq.mapi ( fun i x -> sw.WriteLine(zmw+","+(i+1).ToString()+","+x.Seq.Length.ToString()+x.RQ.ToString()) ) |> ignore

let res = data.CCSReads |> Seq.map outRead |> Seq.length
sw.Close()

let seq1 = data.CCSReads |> Seq.find (fun x -> x.Seq.Count = 136L) 
let seq2 = data.CCSReads |> Seq.find (fun x -> x.Seq.Count = 241L)  
let seq3 = data.CCSReads |> Seq.find (fun x -> x.Seq.Count = 306L)  


let printAln (read:CCSRead) = 
    let alns = data.References |> Seq.map (fun x -> x.AlignSequence(read.Seq)) |> Seq.concat |> Seq.toArray
    Console.WriteLine("Length: " + read.Seq.Count.ToString());
    Console.WriteLine("Seq: "+ read.Seq.ConvertToString());
    alns |> Seq.iter (fun z -> Console.WriteLine(z.ToString())) 


//[<EntryPoint>]
let main args =
[seq1; seq2; seq3] |> List.map printAln |> ignore
    
// Now generate consensus for the third sequence, as peaks don't seem to agree.


let hits = data.CCSReads |> Seq.where (fun x -> x.Seq.Count = 306L) 
                         |> Seq.sortBy (fun z -> -z.SubReads.Count)
                         |> Seq.take 6
                         |> Seq.toArray
hits |> Seq.iter (fun j -> Console.WriteLine(j.SubReads.Count.ToString()))
let mostReads = hits |> Seq.map (fun j -> j.Seq.ConvertToString()) |> Seq.where (fun z -> z.StartsWith("GAAG")) |> Seq.toArray
// Verify top 5 are the same, consider this "correct" 

let allSame i = 
        let char = mostReads.[0].[i]
        mostReads |> Seq.skip 1 |> Seq.exists (fun z -> z.[i] <> char)

let res = {0 .. 305} |> Seq.exists allSame

Console.WriteLine mostReads.[0]



 \|>|>


