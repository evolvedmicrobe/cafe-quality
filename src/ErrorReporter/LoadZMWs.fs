module LoadZMWs

open System
open System.IO
open VariantCaller
open Bio.Util

(* Load the CCS Results from a given directory *)

let ccs_data d_direct = 
    let ccsFiles = (new DirectoryInfo(d_direct)).GetFiles() |> 
                    Seq.where (fun h -> h.Name.EndsWith(".ccs.fastq")) |> 
                    Seq.map (fun u-> u.FullName) |>
                    Seq.toList

    let reference = Path.Combine(d_direct, "../../../References.fna")
    //let reference = "/Users/nigel/git/cafe-quality/data/References.fna"

    // FIXME: This is temporarily hardcoded since I goofed the first runs
    let subReadDirectory = Path.Combine(d_direct, "../../../")
    //let subReadDirectory = "/Users/nigel/CCS_P6_C4/Analysis_Results"
    let subReads = (new DirectoryInfo(subReadDirectory)).GetFiles() |> 
                    Seq.where (fun h -> h.Name.EndsWith(".subreads.fasta.gz")) |> 
                    Seq.map (fun u-> u.FullName) |>
                    Seq.toList



    let qc_exp = new QualityExperiment(null, ccsFiles, subReads, reference)
    Console.WriteLine qc_exp.CCSReads.Count
    qc_exp
