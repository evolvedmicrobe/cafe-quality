namespace PacBio.ConsensusTools

open System.IO

open PacBio.Consensus
open PacBio.IO
open PacBio.FSharp.Utils
open PacBio.FSharp.Utils.ActivePatterns
open PacBio.FSharp.Utils.Logging
open PacBio.FSharp.Utils.SubCommands
open PacBio.Utils
open ManyConsole

type CircularConsensus() as this =
    inherit ConsensusSubCommand ("CircularConsensus", "Generate consensus sequences from single molecules")

    let outputPath = ref "."
    let fofn = ref ""
    let barcodesPath = ref ""
    let minFullPasses = ref 1
    let minPredictedAccuracy = ref 90.0f
    let zmwStart = ref 0
    let zmwCount = ref 0
    let minLength = ref 1
    let maxLength = ref System.Int32.MaxValue
    let doDirectional = ref false
    let doOutputCsv = ref false

    do
        this.optA "o|output=" outputPath "Directory to write results files. Will be created if it doesn't exist."
        this.optA "fofn=" fofn "File of filenames containing bas.h5 files to use."
        this.optA "b|barcodes=" barcodesPath "Barcode .fofn or bc.h5 file. FASTA and FASTQ files will be emitted per-barcode."
        this.optA "minFullPasses=" minFullPasses "Minimum number of complete passes required for CCS reads."
        this.optA "minPredictedAccuracy=" minPredictedAccuracy "Minimum predicted accuracy of CCS reads in percent."
        this.optA "zmwStart=" zmwStart "First ZMW in a range to process."
        this.optA "zmwCount=" zmwCount "Number of ZMW blocks in a range to process."
        this.optA "minLength=" minLength "Minimum length of CCS reads in bases."
        this.optA "maxLength=" maxLength "Maximum length of CCS reads in bases."
        this.optF "directional" (fun _ -> doDirectional := true) "Create separate CCS from forward and reverse subreads. Default = false" 
        this.optF "csv" (fun _ -> doOutputCsv := true) "Output diagnostic information to a CSV file"
#if DEBUG
        this.optF "snrCut=" (fun s -> snrCut := SnrCut.Parse(s)) "Only process ZMWs within a given window of SNR space"
#endif  
        this.AllowsAnyAdditionalArguments () |> ignore

    let paramsFile = "CCSParameters.ini"

    override this.Run(args) =
        this.PrepareOutputDirectory !outputPath

        let processMoviePart basFile ccsH5 fasta fastq csv =
            let basReader = BaseReader.CreateSource basFile
            let pulseFeatures = Set.empty // FIXME: okay for now, but should use PacBio.Data APIs in the future
           
            let range = new ZmwRange(Start = !zmwStart, Count = !zmwCount, Block = 1, Stride = 1)
            let toProcess = basReader.ByHoleNumberRange range

            let stream = new CCSStream ()
           
            try
                let mapFun = (fun x -> (x, stream.Map x))
                let cnt = toProcess |> Seq.parMap mapFun |> Seq.length
                System.Console.WriteLine(cnt)
            with
                | ex ->
                    reraise ()

        let basFiles =
            match (!fofn, args) with
            | Named(path) ->
                if args.Length > 0 then
                    invalidOptf "Unknown arguments: '%s'" <| String.concat " " args
                else
                    Fofn.Filenames path
            | Unnamed(paths) -> paths
            | NonExistent(paths) ->
                invalidOptf "Input file(s) '%s' do(es) not exist -- specify .bas or .bax files on the command line or use the --fofn argument." (String.concat "', '" paths)
            | _ ->
                invalidOpt "Input files unspecified -- specify .bas or .bax files on the command line or use the --fofn argument."

        for basFile in basFiles do
            let path = Path.GetFileName basFile
            let parts = path.Split '.'
            let stem =
                if System.Int32.TryParse (parts.[1], ref 0) then
                    String.concat "." (Seq.take 2 parts)
                else
                    parts.[0]
            let ccsH5 = Path.Combine (!outputPath, stem + ".ccs.h5")
            let fasta = Path.Combine (!outputPath, stem + ".ccs.fasta")
            let fastq = Path.Combine (!outputPath, stem + ".ccs.fastq")
            let csv   = if !doOutputCsv then Path.Combine (!outputPath, stem + ".ccs.csv") else null

            this.logf Info "Processing bas.h5 file: '%s'" basFile

            processMoviePart basFile ccsH5 fasta fastq csv

        int ProcessExitCode.Success
