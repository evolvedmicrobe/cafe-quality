module PacBio.Analysis.Datasets

open System
open PacBio.Align
open PacBio.IO
open PacBio.Utils
open TraceSet

//open PacBio.Common.FSUtils.Utils
//open PacBio.Analysis.CSV

(* Exception type t be used for raising failures that occur
   within Edna and must stop program execution, but are not entirely unexpected
   and therefore should fail gracefully, e.g. if someone tries to train
   without enough data.  The exception should come with a useful error message. *)
exception EdnaError of string


let Log = PacBioLogger.GetLogger()

let logProgress msg =
    Log.Log(LogLevel.INFO, msg)

let logError msg =
    Log.Log(LogLevel.ERROR, msg)
    
let logException (exn : Exception) =
    Log.Log(LogLevel.ERROR, exn.Message)


let zScore (sa : IAlignment) =
    let acc = Math.Max(Math.Min(0.999, sa.Accuracy), 0.2)
    let a = -77.27
    let c = 0.0854
    let sigma = 0.00121
    let L = float sa.TemplateLength
    let d = (a / (L+20.0) + Math.Log(acc / (1.0-acc))  + c) / Math.Sqrt(a*a + 1.0)
    d / sigma
    

type alignmentSource = Internal | SmrtPipe

type ParameterEstimationConfig = {
    Parallel : bool;
    SaveFile : string;
    AnalysisFolder : string option;
    SubSample : int;
    }
    
let baseFind (bmap : char[]) (b : char) = 
    try
        Array.findIndex (fun bc -> bc = b) bmap
    with _ -> failwith(sprintf "Base %O was not found in Basemap: %A" b bmap)

module Seq =
    let has items item =
        Seq.exists (fun i -> i = item) items

let complement c = 
    match c with
        | 'A' -> 'T'
        | 'C' -> 'G'
        | 'G' -> 'C'
        | 'T' -> 'A'
        | _ -> failwith("Invalid base")    
    
    
type Template = { Sequence: char[]; Circular: bool}

(*    
type Trace = 
    abstract member ZmwBases : PacBio.IO.IZmwBases
    abstract member Bases : char[]
    abstract member BaseMap : PacBio.Internal.BaseMap
    abstract member HoleNumber : int
    abstract member Metadata : IMovieMetadata
*)    

type Template with
    member x.IntSequence bmap = x.Sequence |> Array.map (baseFind bmap)
    static member of_string_rc (str:string) (circ:bool) =
        let s = str.ToUpper()
        let tpl = 
            s.ToCharArray() |> Array.filter (Seq.has ['A';'C';'G';'T']) |> Seq.toList |> 
            List.rev |> List.map complement |> List.toArray
        {Sequence = tpl; Circular = circ}
    
    static member to_string_wrap readLength (t : Template) =
        if t.Circular then
            let s = t.Sequence
            let n = int (Math.Max(1.0, Math.Ceiling(float(readLength) / float(s.Length))))
            new String(Array.init (s.Length * n) (fun i -> s.[i%s.Length]))
        else
            new String(t.Sequence)            

(*
// Always return pulses -- handle the different cases of alignments to pulses or bases
let pulsesFromAlignment (trace : Trace) (al : IAlignment) = 
    let bb = trace.Bases
    Array.init al.ReadLength (fun i -> bb.[i + al.ReadStartBase]) |> Array.map (fun b -> b.Pulse)


let baseFilter (n : int) (pulses : #seq<Pulse>) = pulses |> Seq.filter (fun p -> p.Channel = n)
*)

let invBaseMap (tr : Trace) =
    let bmap = tr.Metadata.BaseMap
    let invBmap = bmap |> Seq.mapi (fun channel bse -> (bse, channel)) |> Map.ofSeq
    invBmap


let alignedSequencesIntChunk (al : IAlignment) (tr : Trace) (startCell : AlignCell, endCell : AlignCell) =
    
    let invBmap = invBaseMap tr

    let tpl = 
        let tempString = al.Template.Template       
        let ts = startCell.Template - al.TemplateStartBase
        let tl = endCell.Template - al.TemplateStartBase

        Array.init (tl-ts+1) (fun i -> (int invBmap.[tempString.[ts + i]]) + 1)

            
    let b = (tr.ZmwBases.Base |> Seq.toArray).[startCell.Read .. endCell.Read]
    let rd = b |> Array.map (fun bse -> (int invBmap.[bse]) + 1)
        
    (tpl, rd)


// Get string of the aligned template and read region from an alignment and a trace    
let alignedSequencesIntAll (al : IAlignment) (tr : Trace) =

    let cellRange = al.Cells.[0], al.Cells.[al.Cells.Count-1]
    alignedSequencesIntChunk al tr cellRange

let maxLength = 2500    

let alignedSequencesInt (tr : Trace) (al : IAlignment) =

    let cells = al.Cells |> Seq.toArray
    if cells.Length = 0 then
        ([||],[||])
    else if cells.Length > maxLength then
        let r = new Random(tr.HoleNumber)
        let start = r.Next(cells.Length - maxLength - 1)
        let cellBounds = (cells.[start], cells.[start + maxLength - 1])
        alignedSequencesIntChunk al tr cellBounds
    else
        alignedSequencesIntAll al tr
        

// Movie Uris must always refer to runs, scans, or experiments, and must have a primary analysis report
// tag attached.
type movieSource = Movies of Uri[] | SmrtPipeJob of int


let mkNullable (a:'a) = new Nullable<'a>(a)

let exnKill f =
    try
        Some(f())
    with
        | exn -> Console.WriteLine(exn.Message); Console.WriteLine(exn.StackTrace); None 

