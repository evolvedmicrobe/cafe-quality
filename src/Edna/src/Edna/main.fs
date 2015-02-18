#light
module PacBio.Edna.Main

open System
open System.Linq
open System.Reflection
open System.Diagnostics
open System.IO
open System.Runtime.Serialization
open System.Runtime.Serialization.Formatters.Binary
open System.Text.RegularExpressions
open System.Collections.Generic

open PacBio.Utils
open PacBio.Align

open PacBio.Analysis

open PacBio.Analysis.ErrorEstimation
open PacBio.Analysis.ResultsWriting
open PacBio.Analysis.Datasets

let tester = new PacBio.Hmm.Test.ExpectationMaximization()

open NDesk.Options


let Log = PacBioLogger.GetLogger("Edna")

let logProgress msg =
    Log.Log(LogLevel.INFO, msg)

let logError msg =
    Log.Log(LogLevel.ERROR, msg)
    
let logException (exn : Exception) =
    Log.Log(LogLevel.ERROR, exn.Message)


let fofnFile : string option ref = ref None
let cmpH5Target : string option ref = ref None
let outputFile : string option ref = ref None
let h5out : string option ref = ref None
let maxTraces : int option ref = ref None


// some quick code for finding the cmp.hf file



let _main args =
    let showHelp = ref false

    let popts = 
        [ 
        ("cmp=", "Run a HmmPerTrace workflow for a single cmp.h5 file.", (fun m -> cmpH5Target := Some(m)));
        ("o|output=", "Root filename to write results", (fun o -> outputFile := Some(o)));
        ("h5out=", "Write edna results to an h5 table", (fun o -> h5out := Some(o)));
        ("fofn=", "fofn file", (fun o -> fofnFile := Some(o)));
        ("h|help", "Show help on usage", (fun _ -> showHelp := true));
        ("maxTraces=", "Upper limit on the number of traces that will be processed", (fun arg -> maxTraces := Some(Int32.Parse(arg))));
        ("n|numThreads=", "Maximum number of worker threads to use", (fun arg -> Settings.NumThreads <- Int32.Parse(arg)));
        ]

    try
        let o = new OptionSet()
        popts |> Seq.iter (fun (a,d,f) -> o.Add(a,d, new Action<string>(f)) |> ignore)
        o.Parse <| Environment.GetCommandLineArgs () |> ignore

        Console.WriteLine ("Edna {0}", BuildVersion.Copyright)

        // Process a cmp h5 (in cmph5 order)
        match !cmpH5Target with
            | Some(cmp) -> 
                let cmpFile = if File.Exists(cmp) then cmp else failwith(sprintf "cmp.h5 file does not exist: %s" cmp)
                if String.IsNullOrEmpty(cmpFile) then failwith(sprintf "unrecognized cmp.h5 file: %s" cmp)

                if !fofnFile = None then
                    failwith "Must pass a fofn file name with --fofn argument"

                // Output file
                let outFile =
                    match !outputFile with
                        | None -> Path.Combine(Path.GetDirectoryName(cmpFile), "edna")
                        | Some(f) -> f

                CmpH5Hmm cmpFile outFile (!h5out) (!fofnFile).Value !maxTraces

            | None -> ()

        if (!showHelp) then
            o.WriteOptionDescriptions Console.Out
        0
    with
        | EdnaError msg -> 
            printfn  "%s" ("Cannot run Edna: " + msg)
            1
        | e -> 
            printfn "%s" e.Message
            printfn "%s" e.StackTrace
            1
        
[<STAThread>]
[<EntryPoint>]
let realMain args =
    let r = _main args
    PacBioLogger.Shutdown()
    Environment.Exit(r)
    r