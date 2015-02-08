module PacBio.Analysis.ResultsWriting

open System
open System.Collections.Generic
open System.IO

open PacBio.Align
open PacBio.IO
open PacBio.HDF
open PacBio.Utils

open TraceSet

open PacBio.Hmm.ExpectationMaximization
open PacBio.Hmm.Recursions

open PacBio.Analysis.Datasets
open PacBio.Analysis.CSV


let Log = PacBioLogger.GetLogger()

let logProgress msg =
    Log.Log(LogLevel.INFO, msg)

let logError msg =
    Log.Log(LogLevel.ERROR, msg)
    
let logException (exn : Exception) =
    Log.Log(LogLevel.ERROR, exn.Message)



let setupRunLog outfile =
    let t = Path.GetFullPath(outfile)
    let tRoot = Path.GetFileNameWithoutExtension(t)
    let tp = Path.GetDirectoryName(t)
    let logFile = Path.Combine(tp, tRoot + ".log")
    PacBioLogger.SetupFileLog logFile
    

let setupCentralEndaLog () =
    let logPath = 
        if Environment.OSVersion.Platform <> PlatformID.Win32NT then
            @"/mnt/data3/scratch/Data/edna/log"
        else
            @"\\usmp-data3\scratch\Data\edna\log\"
        
    let now = DateTime.Now
    let fn = Environment.UserName + "-" + (sprintf "%d-%d-%d_%d" now.Month now.Day now.Year now.Second) + "-edna.log"

    let logFile = Path.Combine(logPath, fn)
    PacBioLogger.SetupFileLog logFile

let co i = i :> obj


//let datasetFields (dataset : Dataset) = 
//    dataset.AuxData |> 
//    List.map (fun (name, data) -> (name, data |> co)) |> List.toArray
    

let zmwCode (t : Trace) = t.Metadata.MovieName + "/" + t.HoleNumber.ToString()  
   
let traceFields (t : Trace) = 
    let flds : List<string * obj> = new List<_>();
    let add (fld : string) (f : 'a) = flds.Add(fld,co f)
    
    // Bunch of fields associated with a single trace
    add "Run"           t.Metadata.RunCode
    add "Reports"       t.ReportsFolder
    add "Movie"         (t.Metadata.MovieName)
    add "HoleNumber"    t.HoleNumber
    add "X"             (t.ZmwBases.Zmw.X)
    add "Y"             (t.ZmwBases.Zmw.Y)
    add "Instrument"    (t.Metadata.InstrumentName)
    add "TraceNBases"   (t.ZmwBases.NumBases)
    
    flds


let alignmentFields (t:Trace) (a : IAlignment option) =
    let flds : List<string * obj> = new List<_>();
    let add (fld : string) (f : 'a) = flds.Add(fld,co f)
    
    let ms str = str
    let hasAln (aln : IAlignment) =    
    
        //let n = t.ZmwPulses.StartTime.Count - 1    
        //let startTime = t.ZmwPulses.StartTime.[t.ZmwBases.PulseIndex.[aln.ReadStartBase]]
        //let endTime = t.ZmwPulses.StartTime.[Math.Min(t.ZmwBases.PulseIndex.[aln.ReadStartBase + aln.ReadLength - 1],n)]
        
        let baseWidth = t.ZmwBases.WidthInFrames
        let baseIpd = t.ZmwBases.PreBaseFrames

        let frameTime = float (1.0f / t.Metadata.FrameRate)
        let baseTime = Seq.zip baseWidth baseIpd |> Seq.map (fun (a,b) -> (float (a+b)) * frameTime) |> Seq.toArray

        let baseStartTime = Seq.scan (+) 0.0 baseTime |> Seq.toArray

        let startTime = baseStartTime.[aln.ReadStartBase]
        let endTime = baseStartTime.[aln.ReadStartBase + aln.ReadLength - 1]
        

        let tplName = aln.Template.Name
        add (ms "TemplateName") (if String.IsNullOrEmpty(tplName) then "" else tplName)
        add (ms "TemplateLength") aln.TemplateLength
        add (ms "ReadStart") aln.ReadStartBase
        add (ms "ReadEnd") (aln.ReadStartBase + aln.ReadLength)
        add (ms "TemplateStart") aln.TemplateStartBase
        add (ms "TemplateEnd") (aln.TemplateStartBase + aln.TemplateLength)
        add (ms "StartTime") startTime
        add (ms "StopTime") endTime           
        add (ms "Accuracy") aln.Accuracy
        add (ms "ZScore") (zScore aln)
        add (ms "GlobalRate") ((float aln.TemplateLength) / (endTime - startTime))

        //let p = pulsesFromAlignment t aln
        let alnBaseTimes = baseTime.[aln.ReadStartBase.. aln.ReadStartBase + aln.ReadLength - 1]
        let trimmedBaseTime = Stats.TrimmedMean(alnBaseTimes, 0.05)        
        add (ms "LocalRate") (1.0 / trimmedBaseTime)
    
    let noAln () =
        add (ms "TemplateName") ""
        add (ms "TemplateLength") -1
        add (ms "ReadStart") -1
        add (ms "ReadEnd") -1
        add (ms "TemplateStart") -1
        add (ms "TemplateEnd") -1
        add (ms "StartTime") -1
        add (ms "StopTime") -1        
        add (ms "Accuracy") -1
        add (ms "ZScore") -1  
        add (ms "GlobalRate") -1
        add (ms "LocalRate") -1      
        
    match a with
        | Some(aa) -> if aa <> null && aa.ReadLength > 1 && aa.ReadStartBase >= 0 then hasAln aa else noAln()
        | None -> noAln()
            
    flds


let traceChannelFields (t : Trace, a : IAlignment option, i : int) =
    let flds : List<string * obj> = new List<_>();
    let add (fld : string) (f : 'a) = flds.Add(fld,co f)
    let addset (flds : string list) (f : (Trace*int) -> obj list) = ()
    
    let count f = Array.fold (fun count el -> if f(el) then count + 1 else 0) 0

    //let chp = t.Pulses |> Seq.filter (fun p -> p.Channel = i) |> Seq.toList

    add "Base"              (t.Metadata.BaseMap.[i])
    add "Channel"           (i)
    add "Baseline"          (t.BaselineBias.[i])
    add "Sigma"             (t.BaselineSigma.[i])
    //add "NPulses"           (chp |> Seq.length)    
    
    flds

(*
let totalSignalField (pulses : Pulse[], ch : int) =
    let flds : List<string * obj> = new List<_>();
    let add (fld : string) (f : 'a) = flds.Add(fld,co f)  
    
    add "TotalSignal" (pulses |> Seq.filter (fun p -> p.Channel = ch) |> Seq.map (fun p -> p.Signal) |> Seq.sum)
    flds
*)

let snrSummaryFields (trace: Trace) dtFromChannel ch  = 
    let pkmidToSigma = trace.ZmwBases.Metrics.HQRegionSNR.[ch]    
    let pkmid = pkmidToSigma * float32 trace.BaselineSigma.[ch]

    let totalSignal = (dtFromChannel |> Array.sum) * pkmid * trace.Metadata.FrameRate

    [| 
        "PkmidToSigma", co pkmidToSigma; 
        "Pkmid.Mean",  co pkmid;
        "Pkmid.Median",  co pkmid;
        "TotalSignal", co totalSignal;
        "Pkmid.StdDev", co System.Single.NaN;
        "Pkmid.RStd", co System.Single.NaN;
        "Pkmid.Num", co dtFromChannel.Length
    |]



// Given a function that maps from some datatype to a double array, and a descriptive name,
// build a set of fields that compute stats on the double array
let featureSummaryFields (data : #seq<float32>) (name : string) =    
    let nm stat = name + "." + stat
    let flds = [| nm "Mean"; nm "Median"; nm "StdDev"; nm "RStd"; nm "Num" |]
    
    let stats = new DescriptiveStats(data, true)
    let res = [| co stats.Mean; co stats.Median; co stats.StdDev; co stats.RobustStdDev; co stats.N |]
        
    Seq.zip flds res


let hmmInsertsByChannel(p : ednaParams, obsChannel : int, bmap : char[]) =
    let flds : List<string * obj> = new List<_>();
    let add (fld : string) (f : 'a) = flds.Add(fld,co f)
    
    for cognate = 0 to 3 do
        add (sprintf "InsertOn%c" bmap.[cognate]) p.insert.[obsChannel,cognate]
    done
    
    let nce = ({0..3} |> Seq.map (fun i -> p.insert.[obsChannel,i]) |> Seq.sum) - p.insert.[obsChannel,obsChannel]
    
    add "CognateExtra" p.insert.[obsChannel, obsChannel]
    add "NoncognateExtra" nce
    
    flds

let hmmGenParsByChannel(p : ednaParams, cognateChannel : int) =
    let flds : List<string * obj> = new List<_>();
    let add (fld : string) (f : 'a) = flds.Add(fld,co f)
    
    add "Dark"   p.dark.[cognateChannel]
    add "Merge"  p.merge.[cognateChannel]

    for obsChannel = 0 to 3 do
        add (sprintf "MiscallTo%d" obsChannel) p.miscall.[obsChannel,cognateChannel]
    done
    
    flds


let sliceFeatureByChannel (t:Trace) (al : IAlignment) (feature: 'a array) ch =

    let bases = t.ZmwBases.Base
    let invBM = invBaseMap t
    Seq.init al.ReadLength (fun idx -> al.ReadStartBase + idx) |> Seq.choose (fun idx -> if invBM.[bases.[idx]] = ch then Some(feature.[idx]) else None) |> Seq.toArray



// Write a line for each trace/channel pair
let csvObserver tags filename =
    logProgress (sprintf "Writing trace CSV file: %s" filename)
    let csv = new CSVWriter(filename) :> ITableWriter

    let cleanUp () = csv.Dispose()   

    csv.AddCols(tags |> Seq.map(fun (a,b) -> (a, co b)))
                
    let writeRowInner (region, dists, (tau : float[], missing : float[]), genPars) =
        let (trace,aln) = region
        csv.AddCols(traceFields trace)
        csv.AddCols(alignmentFields trace (Some aln))    

        let brange = [|0 .. trace.Metadata.BaseMap.Length-1|]
    
        try
            //let pulses = pulsesFromAlignment trace aln
            let frameTime = 1.0f / trace.Metadata.FrameRate
            
            let convToTime (t : IList<uint16>) = (t :?> uint16[] |> Array.map (fun v -> (float32 v) * frameTime))

            let dtByCh = sliceFeatureByChannel trace aln (convToTime trace.ZmwBases.WidthInFrames)
            let ipdByCh = sliceFeatureByChannel trace aln (convToTime trace.ZmwBases.PreBaseFrames)
        
            // Write one row per channel
            for i in brange do
                let dtData = (dtByCh i)

                csv.AddCols(hmmInsertsByChannel(genPars,i, trace.Metadata.BaseMap))
                csv.AddCols(hmmGenParsByChannel(genPars,i))
                csv.AddCols(traceChannelFields(trace, Some(aln), i))
                csv.AddCols(snrSummaryFields trace dtData i)
                //csv.AddCols(totalSignalField(pulses, i))
                csv.AddCols(featureSummaryFields dtData "dt")
                csv.AddCols(featureSummaryFields (ipdByCh i) "ipd")
                csv.Row()        
        with e -> 
                logProgress (sprintf "Trace %A failed" trace)
                logProgress (sprintf "Message: %s, Stack Trace: %s" e.Message e.StackTrace)
                

    // Some trace will fail -- don't write a row for these.
    let writeRow v = match v with Some(vv) -> writeRowInner vv | None -> ()

    { new System.IObserver<_> with
        member this.OnCompleted() = cleanUp()
        member this.OnError(e) = cleanUp()
        member this.OnNext(v) = writeRow v }
  
// Write an edna results table to an h5 file
// One row per trace, include just the edna metrics
let cmpTableObserver filename =
    logProgress (sprintf "Writing HMM results to cmp.h5 file: %s" filename)

    // Open h5 file
    let f = HDFFile.Open(filename, FileMode.OpenOrCreate, FileAccess.ReadWrite)
    let ch = new HighLevelChunks(f)
    let ff = ch.File
        
    let width = int64 (ednaParams.starter |> ednaParams.toArray).Length
    let dataSpace = ch.File.CreateDataspace([| 0L; width |], [| -1L; -1L |]);

    // Determine column layout
    let dataType = ch.File.CreateDatatype(typeof<float>)
    ch.File.RemoveChild("Edna") |> ignore
    let dataset = ch.File.CreateDataset("Edna", dataType, dataSpace)

    // Attributes describing the meanings of the columns
    dataset.InsertAttribute("ColumnNames", ednaParams.arrayColumnNames()) |> ignore

    let cleanUp () = f.Dispose()

    let size = ref 0L
    
    let writeArray arr =
        dataset.Extend( [| !size + 1L; width |] )

        let ds = dataset.Dataspace;
        ds.SelectHyperslab([| !size; 0L |], [| 1L; 1L |], [| 1L; width |], null)            
        dataset.Write(arr, ds)

        size := !size + 1L
    
    let writeGoodRow (region, dists, (tau : float[], missing : float[]), genPars) = 
        let arr = ednaParams.toArray genPars
        writeArray arr

    let nans = (Array.create (int width) Double.NaN)
    let writeNans() = writeArray nans

    let writeRow v = match v with Some(vv) -> writeGoodRow(vv) | None -> writeNans()

    { new System.IObserver<_> with
        member this.OnCompleted() = cleanUp()
        member this.OnError(e) = cleanUp()
        member this.OnNext(v) = writeRow v }