﻿//#I @"C:\git\cafe-quality\lib\\"
//#r "PacBio.Utils.dll"
//#r "PacBio.HDF.dll"
//#r "PacBio.IO.dll"
//#r "VariantCaller.dll"
module LoadZMWs

open System;
open System.IO
open VariantCaller
open PacBio.Data
open PacBio.HDF
open PacBio.Utils
open PacBio.IO
open Bio.Util


let loadData = 
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


//[<EntryPoint>]
//let main args = 
//    let data = loadData

//    0


    //let trialFile = direc + "m141008_060349_42194_c100704972550000001823137703241586_s1_p0.1.bax.h5"

    //let data = BaseReader.CreateSource trialFile
    //let zmws = data.ByHoleNumberRange(null) |> Seq.toArray
    //zmws.[0].Metrics.
    //zmws.[0].AdapterHits |> Seq.toArray