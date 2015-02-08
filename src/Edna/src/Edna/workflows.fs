module PacBio.Analysis.ErrorEstimation

open System
open System.Text
open System.Collections.Generic
open System.IO
open System.Diagnostics


open System
open System.Collections.Generic
open System.IO

open PacBio.Align
open PacBio.IO
open PacBio.HDF
open PacBio.Utils

open PacBio.Hmm.ExpectationMaximization
open PacBio.Hmm.Recursions

open PacBio.Analysis.Datasets
open PacBio.Analysis.CSV
open PacBio.Analysis.ResultsWriting

open TraceSet




let Log = PacBioLogger.GetLogger("Edna")

let logProgress msg =
    Log.Log(LogLevel.INFO, msg)

let logError msg =
    Log.Log(LogLevel.ERROR, msg)
    
let logException (exn : Exception) =
    Log.Log(LogLevel.ERROR, exn.Message)


let runOnce src (observers : #seq<IObserver<_>>)  =
    let b = observers |> Seq.toList

    try
        for v in src do
            for o in b do
                o.OnNext v
        
    with e ->
        for i in b do
            i.OnError e

    for i in b do
        i.OnCompleted()

logProgress (sprintf "Will use up to %d threads" Settings.NumThreads)

let innerHmmPerTrace (tags : (string*string) list) (traceSet : TraceSet) (csvOut) (h5Out : string option) =
    
    let getIntAlignment (tr : Trace) =
        try
            let al = tr.FullAlignment
            let alInt = alignedSequencesInt tr al
            Some((tr, al), alInt)
        with
            | _ -> None
            
    let traceRegionsAndAlInt() = traceSet.Traces() |> Seq.map getIntAlignment |> Seq.mapi (fun i v-> (i,v)) 
        
        
    // Fit HMM to data
    logProgress "Estimating 4C HMM Parameters"
    

    let preTrainDataSize = 25
    let bestModelPars =
        // Take the last 25 of [0..max(end, 75)]
        let preTrainData =
            traceRegionsAndAlInt()
            |> Seq.choose snd
            |> Seq.truncate 75
            |> List.ofSeq
            |> List.rev
            |> Seq.truncate preTrainDataSize
            |> Seq.map snd
            |> List.ofSeq
            |> List.rev
        if preTrainData.Length <> preTrainDataSize then  
            let msg = sprintf "There were not enough valid alignments for training. \
                       You needed at least %d alignments.  Please check that the file contains enough." preTrainDataSize
            raise (EdnaError msg)
        estimateTransEmDists transEmDists.starter 0.005 0.03 10 preTrainData     


    let traceProcInner (i, ((tr : Trace, al), intAl))= 
        try
            let sw = new Stopwatch()
            sw.Start()

            let pars = estimateTransEmDists bestModelPars.dists 0.005 0.05 6 [intAl]
            let genPars = transEmDistsToEdnaParams pars
        
            let n = tr.Metadata.BaseMap.Length
            let noMissing = (Array.create n -1.0, Array.create n -1.0)

            let missing = (tr, al) |>  (fun _ -> noMissing)

            sw.Stop()
            if i % 100 = 0 then
                logProgress (sprintf "Trace %d - %d bp - Took %d ms" i (fst intAl).Length (sw.ElapsedMilliseconds))

            Some((tr, al), pars, missing, genPars)
        with
            e -> 
                logProgress (sprintf "Trace %d failed" i)
                logProgress (sprintf "Message: %s, Stack Trace: %s" e.Message e.StackTrace)
                None

    let traceProc (i, v) = match v with Some(vv) -> traceProcInner (i, vv) | None -> None
    let regionAndPars = traceRegionsAndAlInt().ParSelect (fun v -> traceProc v)

    let observers = match h5Out with 
                        | Some(fn) ->  [ (csvObserver tags csvOut); (cmpTableObserver fn) ]
                        | None ->  [ (csvObserver tags csvOut) ]

    runOnce regionAndPars observers
                        

// Find the global 4C HMM model -- put out CSVs containing the HMM parameter estimates, the pulse tagging info,
// and the event counting by template location
let CmpH5Hmm (cmpH5 : string) (outCsv : string) (outH5 : string option) (fofnFile : string) (maxTraces : int Option) =

    // Run logging
    setupRunLog outCsv

    // global logging
    setupCentralEndaLog ()

    //logProgress (sprintf "ExeDir: %s" PacBio.Analysis.Engine.PipelineConfiguration.ExeDir)
    logProgress (sprintf "CmdLine: %s" Environment.CommandLine)
    //logProgress (sprintf "Build: %s" PacBio.Common.Version.VersionUtils.SoftwareVersion)

    logProgress "Starting Hmm CmpH5 analysis"

    logProgress (sprintf "Output file: %s" outCsv)
    
    logProgress (sprintf "Loading CmpH5: %s" cmpH5)

    let traceSet = match maxTraces with
                    | (Some count) -> new TraceSet(cmpH5, fofnFile, count)
                    | None         -> new TraceSet(cmpH5, fofnFile)

    if traceSet.Count > 1000000 then raise (EdnaError "Too many alignments in cmp.h5 file, use --maxTraces")
    logProgress (sprintf "Number of subreads in cmp.h5: %d" traceSet.Count)

    innerHmmPerTrace [] traceSet outCsv outH5
