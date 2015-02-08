module TraceSet

open System
open PacBio.Utils
open PacBio.IO
open PacBio.Hmm.Utils

open System.Collections.Generic

let miniCsvReader f =
    let headers = (System.IO.File.ReadLines f |> Seq.head).Split ','
    let lines = System.IO.File.ReadLines f |> Seq.skip 1

    let headerIndex = headers |> Seq.mapi (fun idx k -> (k.Trim(), idx)) |> Map.ofSeq

    seq {
        for l in lines do
            let chunks = l.Split ','
            yield (fun k -> chunks.[headerIndex.[k]])
        done }


let readBaselineData (bm : Map<char,int>) fofnFile =
    
    //let convertPath = if Environment.OSVersion.Platform = PlatformID.Win32NT then BasCollection.ConvertFofnPath else id
    let convertPath = id

    let resultsDirs = System.IO.File.ReadAllLines(fofnFile) |> Seq.map (convertPath >> System.IO.Path.GetDirectoryName) |> Seq.distinct
    let stsFiles = resultsDirs |> Seq.collect (fun d -> System.IO.Directory.EnumerateFiles(d, "*.sts.csv"))

    let blData = new Dictionary<string, float[]>()
    let sgData = new Dictionary<string, float[]>()

    let readStsCsv f =
        let csvRows = miniCsvReader f

        for r in csvRows do
            let zmwId = (r("Movie")).Trim() + "/" + (r("Zmw")).Trim()

            let bias = Array.zeroCreate 4
            bias.[bm.['T']] <- float (r("BaselineLevel_T"))
            bias.[bm.['G']] <- float (r("BaselineLevel_G"))
            bias.[bm.['A']] <- float (r("BaselineLevel_A"))
            bias.[bm.['C']] <- float (r("BaselineLevel_C"))
            blData.[zmwId] <- bias

            let sigm = Array.zeroCreate 4
            sigm.[bm.['T']] <- float (r("BaselineStd_T"))
            sigm.[bm.['G']] <- float (r("BaselineStd_G"))
            sigm.[bm.['A']] <- float (r("BaselineStd_A"))
            sigm.[bm.['C']] <- float (r("BaselineStd_C"))
            sgData.[zmwId] <- sigm
        done


    stsFiles |> Seq.iter readStsCsv

    (blData, sgData)


type Trace = { reader : TraceSet; ZmwBases : IZmwBases; alignment : IAlnSummary } with
    member x.ReadId = x.ZmwBases.Zmw.Movie.MovieName + "/" + (string x.ZmwBases.Zmw.HoleNumber)
    member x.HoleNumber = x.ZmwBases.Zmw.HoleNumber
    member x.Metadata = x.ZmwBases.Zmw.Movie
    member x.FullAlignment = x.reader.GetAlignment(x.alignment)
    member x.BaselineBias = x.reader.GetBaselineBias(x.ReadId)
    member x.BaselineSigma = x.reader.GetBaselineSigma(x.ReadId)
    member x.ReportsFolder = x.reader.ReportsFolder

and TraceSet(cmpFile, fofn, ?maxTraces0 : int) =
    let cmp = OldCmpH5Reader.CreateReader(cmpFile)
    let bas = BasCollection.FromFofn(fofn)
    let maxTraces = defaultArg maxTraces0 cmp.Alns.Count

    let bm = bas.BaseMap |> Seq.mapi (fun idx bse -> (bse, idx)) |> Map.ofSeq
    let (blData, sigmaData) = readBaselineData bm fofn

    member x.Alns = 
        // random sample, but preserve order
        let rng = new RandomCMWC(42)
        cmp.Alns.TargetSample((fun _ -> true), maxTraces, cmp.Alns.Count, rng) |> Seq.toList

    member x.Count = x.Alns.Length

    member x.Traces() =
        seq {
          for a in x.Alns do
                let zmwBases = bas.GetRead(a)

                let tr = { reader = x; ZmwBases = zmwBases; alignment = a }
                yield tr
            done }

    member x.ReportsFolder = cmp.ReportsFolder
    member x.GetAlignment(aln) = cmp.ReadAlignment aln
    member x.GetBaselineBias id = blData.[id]
    member x.GetBaselineSigma id = sigmaData.[id]