#I  "/Users/nigel/git/cafe-quality/lib/" 
#r "PacBio.Utils.dll"
#r "PacBio.HDF.dll"
#r "PacBio.IO.dll"
#r "VariantCaller.dll"
#r "Bio.dll"
#load "LoadZMWs.fs"

open System;
open System.IO
open System.Collections.Generic
open VariantCaller
open PacBio.HDF
open PacBio.Utils
open Bio.Util
open Bio
open LoadZMWs
open System.Linq

type CCSWriter (fname:string) = 
    let sw = new StreamWriter(fname)
    do sw.WriteLine("ZMW,Reference,NumSubReads,NumErrors,Length,NumIndelErrors,NumSNPErrors")
    member this.Output (read : CCSRead) ( vars: List<Variant>) = 
        if read.AssignedReference <> null && read.SubReads <> null then 
            let indel_cnt = vars |>  Seq.where (fun u  -> u.Type = VariantType.INDEL) |> Seq.length
            let snp_cnt = vars |> Seq.where (fun u -> u.Type = VariantType.SNP) |> Seq.length

            let ref = match read.AssignedReference with
                      | null -> "None"
                      | _ -> read.AssignedReference.RefSeq.ID
            
            let toOut = [| read.ZMWnumber.ToString();
                           read.AssignedReference.RefSeq.ID;
                           read.SubReads.Count.ToString();
                           vars.Count.ToString();
                           read.Seq.Count.ToString();
                           indel_cnt.ToString();
                           snp_cnt.ToString() |]
            let line = String.concat "," toOut
            sw.WriteLine(line)

    member this.Close = sw.Close()

type VariantWriter (fname:string) =     

    let join data =
        data
            |> Seq.map (fun x -> x.ToString())
            |> String.concat ","

    let sw = new StreamWriter(fname)        

    do sw.WriteLine("Ref,Pos,zmw,type,length,num_subreads,homopolymerLength,homopolymerChar,indelSize,indeltype")
        
    member this.Write (read : CCSRead) (variant : Variant )= 

        let vtype = match variant.Type : VariantType with 
                    | VariantType.SNP -> "SNP"
                    | VariantType.INDEL -> "Indel"
                    | VariantType.Complex -> "Complex"
                    | _ -> failwith "not gonna happen"
        
        
        let homopolymerLength = match variant.AtEndOfAlignment with
                                  | true -> "-999"
                                  | false -> match variant with
                                             | :? IndelVariant as indel -> indel.HomopolymerLengthInReference.ToString()
                                             | :? SNPVariant as snp -> "1"
                                             | _ -> failwith "type miss"
        
        let goodIndel = (not variant.AtEndOfAlignment) && variant.Type = VariantType.INDEL
        let indeltype = if goodIndel then
                             System.Enum.GetName(typeof<IndelType>, (variant :?> IndelVariant).InsertionOrDeletion) else
                            "NA"
        let indelLength = if goodIndel then
                            (variant :?> IndelVariant).InsertedOrDeletedBases.Length.ToString() else
                            "NA"

              

        let homoChar = match variant with
                        | :? IndelVariant as indel -> indel.HomopolymerBase.ToString()
                        | :? SNPVariant as snp -> "N"
                        | _ -> failwith "type miss"

        let toWrite : Object[] = [| variant.RefSeq.ID; variant.StartPosition; read.ZMWnumber; vtype; read.Seq.Count; read.SubReads.Count; homopolymerLength; homoChar; indelLength; indeltype|]

        let toOut = join toWrite
        sw.WriteLine(toOut)
    member this.Close = sw.Close()

let vwriter = new VariantWriter("/Users/nigel/git/cafe-quality/data/variants_ratioRef.csv")
let cwriter = new CCSWriter("/Users/nigel/git/cafe-quality/data/ccs_ratioRef.csv")
let mutable totWithVariants = 0;

let vempty = new List<Variant>()

let outputRead (read:CCSRead) =
    if read.AssignedReference <> null then
        let alns = read.AssignedReference.AlignSequence (read.Seq) |> Seq.toArray
        if alns.Length  = 0 then
            cwriter.Output read vempty else
            let best = alns |> Seq.maxBy (fun z -> z.Score)        
            let variants = VariantCaller.CallVariants (best, read.AssignedReference.RefSeq)
            cwriter.Output read variants
            variants |> Seq.iter (fun v -> vwriter.Write read v)





//[<EntryPoint>]
let main args =
    let sw = System.Diagnostics.Stopwatch.StartNew()
    LoadZMWs.ccs_data.CCSReads  |>Seq.iter outputRead
    sw.Stop()
    printfn "%f" sw.Elapsed.TotalMilliseconds
    cwriter.Close
    vwriter.Close

    // Now to ouptu

    Console.WriteLine("Success");
    0
