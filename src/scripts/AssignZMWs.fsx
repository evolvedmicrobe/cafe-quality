#I @"C:\git\cafe-quality\lib\\"
#I  "/Users/nigel/git/cafe-quality/lib/" 
#r "PacBio.Utils.dll"
#r "PacBio.HDF.dll"
#r "PacBio.IO.dll"
#r "VariantCaller.dll"
#r "Bio.dll"
#load "LoadZMWs.fs"

open System;
open System.IO
open VariantCaller
open PacBio.HDF
open PacBio.Utils
open PacBio.IO
open Bio.Util
open Bio

let data = LoadZMWs.ccs_data



let alnLength = data.References |> Seq.map (fun ref -> ref.RefSeq.Count) |> Seq.toArray

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
    let alns = data.References |> Seq.map (fun x -> x.AlignSequence(read.Seq)) |> Seq.toArray
    alns |> Seq.map (fun z -> Console.Write(z.ToString()))

[seq1; seq2; seq3] |> List.map printAln |> ignore