module LoadZMWs

open System
open System.IO
open VariantCaller
open Bio.Util


let ccs_data = 
    let runninMono = Bio.CrossPlatform.Environment.RunningInMono
    let d_direct = if runninMono  then 
                        "/Users/nigel/git/cafe-quality/data/" else
                        @"C:\git\cafe-quality\data\"

    let direc = if runninMono then
                    "/Users/nigel/CCS_P6_C4/Analysis_Results/"
                else
                    @"C:\CCS_P6_C4\Analysis_Results\\"
                   
    let baxFiles = (new DirectoryInfo(direc)).GetFiles() |>
                    Seq.where   (fun z -> z.Name.EndsWith(".bax.h5")) |>
                    Seq.map (fun u -> u.FullName) |>
                    Seq.toList

    let ccsFiles = (new DirectoryInfo(d_direct)).GetFiles() |> 
                    Seq.where (fun h -> h.Name.EndsWith(".ccs.fasta.gz")) |> 
                    Seq.map (fun u-> u.FullName) |>
                    Seq.toList

    let reference = Path.Combine(d_direct, "References.fna")

    let subReads = (new DirectoryInfo(d_direct)).GetFiles() |> 
                    Seq.where (fun h -> h.Name.EndsWith(".subreads.fasta.gz")) |> 
                    Seq.map (fun u-> u.FullName) |>
                    Seq.toList



    let qc_exp = new QualityExperiment(baxFiles, ccsFiles, subReads)
    Console.WriteLine qc_exp.CCSReads.Count
    qc_exp

 

//[<EntryPoint>]
//let main args = 
//    let data = loadData
//    
//    0


