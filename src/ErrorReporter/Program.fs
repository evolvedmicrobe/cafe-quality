open System
open System.IO
open System.Collections.Generic
open VariantCaller
open Bio.Util
open Bio
open LoadZMWs
open System.Linq

type CCSWriter (fname:string) = 
    let sw = new StreamWriter(fname)
    do sw.WriteLine("Movie,ZMW,Reference,NumSubReads,NumErrors,Length,NumIndelErrors,NumSNPErrors")
    member this.Output (read : CCSRead) ( vars: List<Variant>) = 
        if read.AssignedReference <> null && read.SubReads <> null then 
            let indel_cnt = vars |>  Seq.where (fun u  -> u.Type = VariantType.INDEL) |> Seq.length
            let snp_cnt = vars |> Seq.where (fun u -> u.Type = VariantType.SNP) |> Seq.length

            let ref = match read.AssignedReference with
                      | null -> "None"
                      | _ -> read.AssignedReference.RefSeq.ID
            
            let toOut = [| read.Movie;
                           read.ZMWnumber.ToString();
                           read.AssignedReference.RefSeq.ID;
                           read.SubReads.Count.ToString();
                           vars.Count.ToString();
                           read.Seq.Count.ToString();
                           indel_cnt.ToString();
                           snp_cnt.ToString() |]
            let line = String.concat "," toOut
            sw.WriteLine(line)
            sw.Flush()

    member this.Close = sw.Close()

type VariantWriter (fname:string) =     

    let join data =
        data
            |> Seq.map (fun x -> x.ToString())
            |> String.concat ","

    let sw = new StreamWriter(fname)        

    do sw.WriteLine("Ref,Pos,zmw,type,length,homopolymerLength,homopolymerChar,indelSize,indeltype")
        
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

        let toWrite : Object[] = [| variant.RefSeq.ID; variant.StartPosition; read.ZMWnumber; vtype; read.Seq.Count; homopolymerLength; homoChar; indelLength; indeltype|]

        let toOut = join toWrite
        sw.WriteLine(toOut)
        sw.Flush()
    member this.Close = sw.Close()


let mutable totWithVariants = 0;

let vempty = new List<Variant>()
let mutable outCnt = 0
let outputRead (vwriter:VariantWriter) (cwriter:CCSWriter) (read:CCSRead) =
    if read.AssignedReference <> null then
        outCnt <- outCnt + 1
        if outCnt % 50 = 0 then Console.WriteLine ("Outputting read number: " + outCnt.ToString())
        let alns = read.AssignedReference.AlignSequence (read.Seq) |> Seq.toArray
        if alns.Length  = 0 then
            cwriter.Output read vempty else
            let best = alns |> Seq.maxBy (fun z -> z.Score)        
            let variants = VariantCaller.CallVariants (best, read.AssignedReference.RefSeq)
            cwriter.Output read variants
            variants |> Seq.iter (fun v -> vwriter.Write read v)


[<EntryPoint>]
let main args =
    let sw = System.Diagnostics.Stopwatch.StartNew()
    let direc = args.[0]
    let prefix = args.[1]
    let fname_prefix = Path.Combine(direc,prefix)
    let vwriter = new VariantWriter( fname_prefix + "_all_variants.csv")
    let cwriter = new CCSWriter(fname_prefix + "_read_report.csv")
    let data = LoadZMWs.ccs_data direc
    let outputter = outputRead vwriter cwriter
    data.CCSReads  |> Seq.iter outputter
    sw.Stop()
    printfn "%f" sw.Elapsed.TotalMilliseconds
    cwriter.Close
    vwriter.Close
    Console.WriteLine("Success");
    0
