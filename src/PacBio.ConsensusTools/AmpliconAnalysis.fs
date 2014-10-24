
namespace PacBio.ConsensusTools

open System.IO
open System.Collections.Generic

open PacBio.Consensus
open PacBio.IO
open PacBio.IO.Fasta
open PacBio.FSharp.Utils
open PacBio.FSharp.Utils.ActivePatterns
open PacBio.FSharp.Utils.Logging
open PacBio.FSharp.Utils.SubCommands
open PacBio.Utils

type AmpliconAnalysis() as this =
    inherit ConsensusSubCommand ("AmpliconAnalysis", "Assemble and call consensus for a collection of diploid amplicons")

    let minPredictedAccuracy = ref 0.95f
    let fofn = ref ""
    let outputPath = ref "."
    let barcodesPath = ref ""
    let barcodes = ref ([] : string list)
    let whiteList = ref ""
    let maxReads = ref 2000
    let maxPhasingReads = ref 500
    let minLength = ref 3000
    let maxLength = ref 0
    let minReadScore = ref 0.75f
    let minSnr = ref 3.0f
    let prioritySample = ref false
    let maxClusters = ref 0
    let doClustering = ref true
    let doPhasing = ref true
    let minSplitScore = ref 6.0
    let minSplitFraction = ref 0.1
    let minSplitReads = ref 20
    let minAvgBarcodeScore = ref 0
    let minTotalBarcodeScore = ref 0
    let sampleName = ref "0"
    let clusterInflation = ref 2.0f
    let doChimeraFilter = ref true
    let chimeraScoreThreshold = ref 1.0f
    let ignoreEnds = ref 0
    let trimEnds = ref 0
    let takeN = ref 0
    let doFastxByBarcode = ref false
    let doSubreadDebug = ref false
    let doConvergenceFilter = ref false
    let doCCS = ref false
    let forcePar = ref false
    let precision = ref 5

    do
        // set a bunch of the variables for CCS mode
        let enableCCS () =
            doCCS := true
            doClustering := false
            doPhasing := false
            doChimeraFilter := false
            minLength := 0
            maxLength := 0
            minSnr := 0.0f
            minReadScore := 0.0f
            minSplitReads := 1

        this.optF "circularConsensus" enableCCS "Build circular consensus sequences from subreads per ZMW. Default = false"
        this.optA "o|output=" outputPath "Directory to write results files. Will be created if it doesn't exist."
        this.optA "f|fofn=" fofn "File of filenames containing bas.h5 files to use."
        this.optA "b|barcodes=" barcodesPath "Barcode .fofn. AmpliconAnalysis will work on each barcode independently."
        this.optF "doBc=" (fun v -> barcodes := v::!barcodes) "Use this flag to specify a subset of barcodes to process."
        this.optA "whiteList=" whiteList "A list of ReadIds to use, in FASTA format read headers or as text, one per line. An example ReadId is '<movie_name>/<hole_number>'."
        this.optA "r|maxReads=" maxReads "Number of subreads to use."
        this.optA "maxPhasingReads=" maxPhasingReads "Number of subreads to use for phasing/consensus of each coarse cluster."
        this.optA "l|minLength=" minLength "Minimum subread length to use."
        this.optA "L|maxLength=" maxLength "Maximum subread length to use. Supply a value <1 to disable."
        this.optA "s|minReadScore=" minReadScore "Minimum read score to use."
        this.optA "minPredictedAccuracy=" minPredictedAccuracy "Minimum predicted accuracy below which a haplotype is considered 'noise'."
        this.optA "minSnr=" minSnr "Minimum SNR of reads."
        this.optF "prioritySample" (fun () -> prioritySample := true) "Priority sample reads rather than reservoir sampling them. Default = false"
        this.optA "maxClusters=" maxClusters "Return the sequences of the best clusters, up to this number. Set to 0 for unlimited."
        this.optF "noClustering" (fun () -> doClustering := false) "Deactivate the clustering step. Read go directly into phasing and consensus. Default = false"
        this.optF "noPhasing" (fun () -> doPhasing := false) "Deactivate the phasing step. Call a single consensus sequence for each read cluster. Default = false"
        this.optA "minSplitScore=" minSplitScore "Global quiver score improvement required to split a haplotype. Set to a value > 10000 to disable phasing entirely."
        this.optA "minSplitFraction=" minSplitFraction "Minimum fraction of reads favoring the minor phase required to split a haplotype."
        this.optA "minSplitReads=" minSplitReads "Minimum number reads favoring the minor phase required to split a haplotype."
        this.optA "minBarcodeScore=" minAvgBarcodeScore "Minimum average barcode score to require of barcoded subreads."
        this.optA "minTotalBarcodeScore=" minTotalBarcodeScore "Minimum total barcode score to require of barcoded subreads."
        this.optA "sampleName=" sampleName "Descriptive name to add to FASTA header. Only used in non-barcoded mode."
        this.optA "clusterInflation=" clusterInflation "Markov clustering inflation parameter, larger values give more clusters."
        this.optF "noChimeraFilter" (fun () -> doChimeraFilter := false) "Deactivate the chimera filter. All haplotypes will be returned as normal sequences. Default = false"
        this.optA "chimeraScoreThreshold=" chimeraScoreThreshold "Deactivate the chimera filter. All haplotype will be returned as normal sequences."
        this.optA "i|ignoreEnds=" ignoreEnds "When splitting, ignore N bases at the ends. Use to prevent excessive splitting caused by degenerate primers."
        this.optA "t|trimEnds=" trimEnds "Trim N bases at each end of each consensus sequence. Used to remove possible errors induced by degenerate primers."
        this.optF "x|fastxByBarcode" (fun () -> doFastxByBarcode := true) "Write consensus sequences out to separate files by barcode. Default = false"
        this.optA "takeN=" takeN "Report only the top N sequences (by number of supporting subreads) from each coarse cluster, the remainder will be considered 'noise'. Supply a value <1 to disable."
        this.optF "subreadDebug" (fun () -> doSubreadDebug := true) "Write a csv file containing subread summary info (including barcode data if applicable), for debugging purposes. Default = false"
        this.optF "convergenceFilter" (fun () -> doConvergenceFilter := true) "Enable the the convergence filter. Only haplotypes whose Quiver consensus converges will be returned as normal sequences. Default = false"
        this.optF "useTheForce" (fun () -> forcePar := true) "Force threaded parallelism across barcodes. (WARNING: MAY CRASH WITH OUT OF MEMORY ERROR)"

        this.AllowsAnyAdditionalArguments () |> ignore

    override this.Run(args) =
        MultiTemplateConsensus.minSplitFraction <- !minSplitFraction
        MultiTemplateConsensus.splitThreshold <- !minSplitScore
        MultiTemplateConsensus.minSplitReads <- !minSplitReads
        MultiTemplateConsensus.ignoreEnds <- !ignoreEnds
        FineClustering.DoPhasing <- !doPhasing

        this.PrepareOutputDirectory !outputPath

        let basCollection =
            match (!fofn, args) with
            | Named(path) ->
                if args.Length > 0 then
                    invalidOptf "Unknown arguments: '%s'" <| String.concat " " args
                else
                    BasCollection.FromFofn path
            | Unnamed(paths) -> BasCollection paths
            | NonExistent(paths) ->
                invalidOptf "Input file(s) '%s' do(es) not exist -- specify .bas or .bax files on the command line or use the --fofn argument." (String.concat "', '" paths)
            | _ ->
                invalidOpt "Input files unspecified -- specify .bas or .bax files on the command line or use the --fofn argument."

        // try to load the chemistries to force an exception here rather than later
        let parametersFile = if !doCCS then "CCSParameters.ini" else "QuiverParameters.ini"
        let pulseFeatures = Set.empty // FIXME: okay now, fix when moving to PacBio.Data APIs
        let chemistries = try set basCollection.Chemistries with | _ -> Set.empty
        use scorerConfig = this.LoadQuiverConfig parametersFile pulseFeatures chemistries

        // if doSubreadDebug then write out some properties about all subreads in the collection
        let writeReadsFile (bas : BasCollection) = do
            // TODO: convert this to use F# data providers
            let path = Path.Combine (!outputPath, "amplicon_assembly_subread_debug.csv")
            use writer = new CsvWriter (path)
            for sr in bas.Subreads do
                writer.AddCol ("SubreadID", sr.SubreadId)
                writer.AddCol ("Length", sr.Region.Length)
                writer.AddCol ("AdapterBefore", sr.Region.AdapterHitBefore)
                writer.AddCol ("AdapterAfer", sr.Region.AdapterHitAfter)
                writer.AddCol ("ReadScore", sr.ReadScore)
                writer.AddCol ("MinSNR", sr.MinSnr)
                writer.Row ()

        if !doSubreadDebug then
            writeReadsFile basCollection

        // do the actual amplicon assembly here
        let doAmpliconAnalysis (allReads : seq<Subread>) (name : string) =
            if name <> !sampleName then do
                this.logf Info "Barcode '%s', filtering from %d total subreads..." name (Seq.length allReads)
            else
                this.logf Info "Filtering from %d total subreads..." (Seq.length allReads)
            
            let readFilter (sr : Subread) =
                sr.Region.Length >= !minLength &&
                (!maxLength < 1 || sr.Region.Length <= !maxLength) &&
                sr.ReadScore >= !minReadScore &&
                sr.MinSnr >= !minSnr
            let filtered = allReads |> Seq.filter readFilter

            if name <> !sampleName then do
                this.logf Info "Barcode '%s', sampling from %d post-filter subreads..." name (Seq.length filtered)
            else
                this.logf Info "Sampling from %d post-filter subreads..." (Seq.length filtered)

            let reads =
                if !prioritySample then
                    filtered
                    |> Seq.toList
                    |> List.sortWith (fun a b -> -a.ReadScore.CompareTo(b.ReadScore))
                    |> Seq.truncate !maxReads
                    |> Seq.toArray
                else
                    filtered.ReservoirSample !maxReads
                    |> Seq.toList
                    |> List.sortWith (fun a b -> a.SubreadId.CompareTo(b.SubreadId))
                    |> Seq.toArray

            if name <> !sampleName then do
                this.logf Info "Barcode '%s', coarse clustering %d subreads..." name reads.Length
            else
                this.logf Info "Coarse clustering %d subreads..." reads.Length

            let clusters =
                CoarseClustering.ClusterReads (reads, !clusterInflation, null, !doClustering)
                |> Seq.filter (fun c -> c.Length >= !minSplitReads)
                |> Seq.sortBy (fun c -> -c.Length)
                |> if !maxClusters > 0 then Seq.truncate !maxClusters else id
                |> Seq.mapi (fun i c -> (i, c))
                |> Seq.toArray

            let numClusteredReads = Seq.sumBy (fun (_, s : Subread[]) -> s.Length) clusters
            if name <> !sampleName then do
                this.logf Info "Barcode '%s', fine clustering %d subreads..." name numClusteredReads
            else
                this.logf Info "Fine clustering %d subreads..." numClusteredReads

            let phaseCluster (i, c) =
                FineClustering.FineClusteringAndConsensus (i, c, scorerConfig, name, !maxPhasingReads, !takeN > 0)

            let results = clusters |> Seq.map phaseCluster |> Seq.concat |> Seq.toArray

            // if we're not phasing and not clustering, then don't bother with the chimera filter
            if !doChimeraFilter && (!doPhasing || !doClustering) then
                let chimeraInput = results |> Seq.map (fun r -> (r.Coverage, r.FastaName, r.Sequence))
                let chimeraResults = ChimeraDetector.findChimeras chimeraInput !chimeraScoreThreshold |> snd
                let resDict = dict [ for r in results do yield (r.FastaName, r) ]
                for chimera in chimeraResults do
                    resDict.[chimera.name].ChimeraResult <- chimera
                resDict.Values |> Seq.sortBy (fun r -> -r.Coverage) |> Seq.toArray
            else
                results
        // End of doAmpliconAnalysis enclosure

        // output functions
        let fastaName (cr : ConsensusResult) =
            if !doCCS then cr.BarcodeName else cr.FastaName

        let writeFastas prefix (results : #seq<ConsensusResult>) = do
            let fastaPath = Path.Combine (!outputPath, prefix + ".fasta")

            File.Delete fastaPath

            use fastaWriter = new SimpleFASTAWriter (fastaPath)
            for cr in results do
                if !trimEnds > 0 then
                    let trim = !trimEnds
                    if cr.Sequence.Length > (2 * trim) then
                        fastaWriter.WriteEntry (fastaName cr, cr.Sequence.[trim..cr.Sequence.Length-trim-1])
                else
                    fastaWriter.WriteEntry (fastaName cr, cr.Sequence)
                    
        let writeFastqs prefix (results : #seq<ConsensusResult>) = do
            let fastqPath = Path.Combine (!outputPath, prefix + ".fastq")

            File.Delete fastqPath

            let toQVs = Seq.map (fun (s : ConsensusBaseStats) -> uint32 s.ConsensusQV) >> Seq.toArray
            use fastqWriter = new SimpleFASTQWriter (fastqPath)
            for cr in results do
                let entry =
                    if !trimEnds > 0 then
                        let trim = !trimEnds
                        if cr.Sequence.Length <= 2 * trim then
                            None
                        else
                            let sequence = cr.Sequence.[trim..cr.Sequence.Length-trim-1]
                            let qvs = cr.Stats.[trim..cr.Stats.Length-trim-1] |> toQVs
                            Some <| FASTQEntry (fastaName cr, sequence, qvs)
                    else
                        Some <| FASTQEntry (fastaName cr, cr.Sequence, cr.Stats |> toQVs)
                match entry with
                | Some(e) -> fastqWriter.WriteEntry e
                | None -> ()

        let writeNonBarcodeSequences prefix (results : #seq<ConsensusResult>) = do
            writeFastas prefix results
            writeFastqs prefix results

        let writeBarcodeSequences bcNames prefix (results : #seq<ConsensusResult>) = do
            let fastaDirectoryName = prefix + "_fasta"
            let fastqDirectoryName = prefix + "_fastq"
            let fastaDirectory = Path.Combine (!outputPath, fastaDirectoryName)
            let fastqDirectory = Path.Combine (!outputPath, fastqDirectoryName)
            this.PrepareOutputDirectory fastaDirectory
            this.PrepareOutputDirectory fastqDirectory
            for bc in bcNames do
                let fastaPrefix = Path.Combine (fastaDirectoryName, bc)
                let fastqPrefix = Path.Combine (fastqDirectoryName, bc)
                let bcResults = results |> Seq.filter (fun cr -> cr.BarcodeName = bc)
                writeFastas fastaPrefix bcResults
                writeFastqs fastqPrefix bcResults
        // End of output functions

        if !doCCS then
            let isFullPass (subread : Subread) =
                subread.Region.AdapterHitBefore = true && subread.Region.AdapterHitAfter = true

            let doCircularConsensus movie holeNumber =
                let subreads = basCollection.GetSubreadsForZmw (movie, holeNumber)
                let numPasses = subreads |> Seq.filter isFullPass |> Seq.length
                doAmpliconAnalysis subreads (sprintf "%s/%d/ccs/%d" movie holeNumber numPasses)

            let isGoodResult (cr : ConsensusResult) =
                let isAcc = cr.PredictedAccuracy >= !minPredictedAccuracy
                let isConv = (not !doConvergenceFilter) || cr.ConsensusConverged
                if not isAcc then
                    this.logf Info "Omitting '%s', insufficiently accurate: %.3f" cr.BarcodeName cr.PredictedAccuracy
                isAcc && isConv

            let results =
                basCollection.Movies
                |> Seq.map (fun m -> basCollection.HoleNumber m |> Seq.parMap (fun hn -> doCircularConsensus m hn))
                |> Seq.concat
                |> Seq.concat
                |> Seq.filter isGoodResult

            let writeFastxs prefix (results : #seq<ConsensusResult>) =
                let fastaPath = Path.Combine (!outputPath, prefix + ".fasta")
                let fastqPath = Path.Combine (!outputPath, prefix + ".fastq")

                File.Delete fastaPath
                File.Delete fastaPath

                let toQVs = Seq.map (fun (s : ConsensusBaseStats) -> uint32 s.ConsensusQV) >> Seq.toArray
                use fastaWriter = new SimpleFASTAWriter (fastaPath)
                use fastqWriter = new SimpleFASTQWriter (fastqPath)

                for cr in results do
                    fastaWriter.WriteEntry (fastaName cr, cr.Sequence)
                    fastqWriter.WriteEntry <| FASTQEntry (fastaName cr, cr.Sequence, cr.Stats |> toQVs)

            writeFastxs "circular_consensus" results

        else
            // return the consensus results and how to write out the sequences
            let results, writeSequences =
                match !barcodesPath with
                // do barcoded analysis
                | NonEmpty(bcPath) ->
                    let reader = new BarcodeReader (bcPath)
                    let bcNames =
                        if List.isEmpty !barcodes then
                            reader.BarcodeNames
                        else
                            let valid = Set.ofArray reader.BarcodeNames
                            !barcodes |> Seq.filter (fun n -> Set.contains n valid) |> Seq.toArray

                    let resultsForBarcode bc =
                        let reads =
                            reader.Calls
                            |> Seq.filter (fun c -> c.BarcodeName = bc)
                            |> Seq.filter (fun c -> c.AvgScore >= !minAvgBarcodeScore)
                            |> Seq.filter (fun c -> c.Score >= !minTotalBarcodeScore)
                            |> Seq.map (fun c -> basCollection.GetSubreadsForZmw(c.ReadId))
                            |> Seq.concat
                        this.logf Info "Barcode '%s', starting with %d barcoded subreads..." bc (Seq.length reads) 
                        doAmpliconAnalysis reads bc

                    // disable parallelism if phasing -- we might run out of memory
                    let mapBarcode = if !forcePar || not !doPhasing then Seq.parMap else Seq.map

                    let results = bcNames |> mapBarcode resultsForBarcode |> Seq.concat |> Seq.toArray

                    let writeSequences =
                        if !doFastxByBarcode then
                            writeBarcodeSequences bcNames
                        else
                            writeNonBarcodeSequences

                    (results, writeSequences)

                // do non-barcoded analysis (possibly whitelisted)
                | _ ->
                    if not <| List.isEmpty !barcodes then
                        invalidOpt "Barcodes cannot be used without also specifying a .fofn of barcode h5 files."

                    let reads =
                        match !whiteList with
                        | NonEmpty(wlPath) ->
                            this.logf Info "Using file '%s' as a ZMW whitelist" wlPath
                            let fastaExts = [|".fasta"; ".fsta"; ".fa"; ".fas"|]
                            let fastqExts = [|".fastq"; ".fstq"; ".fq"; ".faq"|]
                            let wlExt = Path.GetExtension wlPath
                            let isExt = fun e -> e = wlExt
                            let readIds =
                                if Array.exists isExt fastaExts then
                                    use reader = new Fasta.SimpleFASTAReader (wlPath)
                                    [| for e in reader -> Subreads.getSubreadZmw e.Header |]
                                elif Array.exists isExt fastqExts then
                                    this.log Error "Whitelist from FASTQ not implemented - bug lhepler@pacificbiosciences.com."
                                    new System.NotImplementedException () |> raise
                                    [||]
                                else
                                    [| for e in File.ReadLines wlPath -> Subreads.getZmw e |]
                            let uniqueReadIds = Seq.distinct readIds |> Seq.toArray
                            this.logf Info "Found %d unique whitelisted ZMWs" uniqueReadIds.Length
                            seq { for readId in Set uniqueReadIds do
                                    for read in basCollection.GetSubreadsForZmw readId do yield read }
                        | _ ->
                            this.log Info "No whitelist file specified, using all available subreads" 
                            basCollection.Subreads

                    let results = doAmpliconAnalysis reads !sampleName

                    (results, writeNonBarcodeSequences)
            // End of results enclosure

            // write the results to the various files (usually the last step)
            let writeResults (results : ConsensusResult[]) = do
                this.log Info "Writing results..."

                let isGoodResult (cr : ConsensusResult) =
                    (String.empty cr.DuplicateOf) &&
                    cr.ChimeraResult.chimeraScore < !chimeraScoreThreshold &&
                    cr.PredictedAccuracy >= !minPredictedAccuracy &&
                    // if the convergence filter is off or the sequence converged
                    ((not !doConvergenceFilter) || cr.ConsensusConverged) &&
                    // if the "takeN" filter is off or the sequence is in the top N
                    ((!takeN <= 0) || cr.PhaseIndex < !takeN)

                this.logf Info "Writing %i consensus sequences to FASTA and FASTQ..." results.Length
                results |> Seq.filter isGoodResult |> writeSequences "amplicon_analysis"
                results |> Seq.filter (isGoodResult >> not) |> writeSequences "amplicon_analysis_chimeras_noise"

                let writePerBaseCsvResults prefix (results : #seq<ConsensusResult>) =
                    this.log Info "Writing per-base results to CSV..."
                    let path = Path.Combine (!outputPath, prefix + ".csv")
                    use writer = new CsvWriter (path)
                    for cr in results do
                        writer.AddCol ("BarcodeName", cr.BarcodeName)
                        writer.AddCol ("FastaName", cr.FastaName)
                        writer.AddCol ("CoarseCluster", cr.CoarseClusterName)
                        writer.AddCol ("Phase", cr.FineClusterName)
                        writer.AddCol ("TotalCoverage", cr.Coverage)
                        writer.AddCol ("SequenceLength", cr.Stats.Length)

                        for i = 0 to cr.Stats.Length - 1 do
                            writer.AddCol ("Position", i)
                            writer.AddCol ("Base", cr.Stats.[i].Base)
                            writer.AddCol ("QV", cr.Stats.[i].ConsensusQV)
                            writer.AddCol ("Coverage", cr.Stats.[i].Coverage)
                            writer.Row ()

                let writePerReadCsvResults prefix (results : #seq<ConsensusResult>) =
                    this.log Info "Writing per-read summary information to CSV..."
                    let path = Path.Combine (!outputPath, prefix + ".csv")
                    use writer = new CsvWriter (path)
                    for cr in results do
                        writer.AddCol ("BarcodeName", cr.BarcodeName)
                        writer.AddCol ("FastaName", cr.FastaName)
                        writer.AddCol ("CoarseCluster", cr.CoarseClusterName)
                        writer.AddCol ("Phase", cr.FineClusterName)
                        writer.AddCol ("TotalCoverage", cr.Coverage)
                        writer.AddCol ("SequenceLength", cr.Stats.Length)
                        writer.AddCol ("PredictedAccuracy", cr.PredictedAccuracy)
                        writer.AddCol ("ConsensusConverged", cr.ConsensusConverged)
                        writer.AddCol ("DuplicateOf", cr.DuplicateOf)
                        writer.AddCol ("NoiseSequence", isGoodResult cr |> not)

                        let isChimera = cr.ChimeraResult.chimeraScore >= !chimeraScoreThreshold

                        writer.AddCol ("IsChimera", isChimera)
                        writer.AddCol ("ChimeraScore", cr.ChimeraResult.chimeraScore)
                        writer.AddCol ("ParentSequenceA", if isChimera then cr.ChimeraResult.sequenceA else "")
                        writer.AddCol ("ParentSequenceB", if isChimera then cr.ChimeraResult.sequenceB else "")
                        writer.AddCol ("CrossoverPosition", if isChimera then cr.ChimeraResult.crossoverPosition else -1)
                        writer.Row ()

                let writeSubreadMappings prefix (results : ConsensusResult[]) =
                    this.log Info "Writing subread mapping scores to CSV..."
                    let path = Path.Combine (!outputPath, prefix + ".csv")
                    use writer = new CsvWriter (path)
                    let subreadIds = 
                        let allSubreads = seq {
                            for cr in results do 
                                for id in cr.SubreadIds do 
                                    yield id
                            }
                        allSubreads 
                        |> Seq.distinct 
                        |> Seq.toList 
                        |> List.sortWith Comparer.subreadIds
                        |> List.toArray
                    for id in subreadIds do
                        writer.AddCol ("SubreadId", id)
                        for cr in results do
                            let score =
                                try
                                    let index = Array.findIndex ((=) id) cr.SubreadIds
                                    let rawScore = cr.MappingScores.[index]
                                    System.Math.Round (float rawScore, 5)
                                with
                                | :? KeyNotFoundException -> 0.
                                | _ as ex -> raise ex
                            writer.AddCol (cr.FastaName, score)
                        writer.Row ()

                let writeMappings prefix (results : ConsensusResult[]) =
                    // Declare the function we will use to write mapping data to file
                    let writeMappingCsv filename columnName (weightTuples : #seq<string * 'a[]>) =
                        let path = Path.Combine (!outputPath, filename)
                        use writer = new CsvWriter (path)
                        for (sequenceId, weights) in weightTuples do
                            writer.AddCol (columnName, sequenceId)
                            for (i, cr) in Seq.enumerate results do
                                writer.AddCol (cr.FastaName, weights.[i])
                            writer.Row ()

                    // Build an array of mapping scores for each subread
                    this.log Info "Calculating per-subread mapping scores..."
                    let subreadIds = 
                        let allSubreads = seq {
                            for cr in results do 
                                for id in cr.SubreadIds do 
                                    yield id
                            }
                        allSubreads 
                        |> Seq.distinct 
                        |> Seq.toList 
                        |> List.sortWith Comparer.subreadIds
                        |> List.toArray
                    let subreadWeights =
                        let consensusWeights id =
                            let weight (cr : ConsensusResult) =
                                try
                                    let index = Array.findIndex ((=) id) cr.SubreadIds
                                    let rawScore = cr.MappingScores.[index]
                                    System.Math.Round (float rawScore, !precision)
                                with
                                | :? KeyNotFoundException -> 0.
                                | _ as ex -> raise ex
                            results |> Array.map weight
                        subreadIds |> Array.map (fun id -> (id, consensusWeights id))

                    // Write the subread mapping scores to file as-is
                    this.log Info "Writing per-subread mapping scores to CSV..."
                    writeMappingCsv (prefix + "_subreads.csv") "SubreadId" subreadWeights

                    // Identify all unique ZMWs and average the scores of their subreads
                    this.log Info "Calculating per-Zmw mapping scores..."
                    let zmwWeights =
                        let subreadsByZmw = 
                            subreadIds 
                            |> Seq.groupBy Subreads.getSubreadZmw
                            |> Map.ofSeq
                        let subreadWeightsMap = subreadWeights |> Map.ofSeq
                        let inline (++) a b = Array.map2 (+) a b
                        let zeroArray = Array.create results.Length 0.0
                        let zmwWeight zmw =
                            let subreads = Map.find zmw subreadsByZmw
                            let count = float (Seq.length subreads)
                            subreads 
                            |> Seq.map (fun s -> Map.find s subreadWeightsMap)
                            |> Seq.fold (++) zeroArray
                            |> Seq.divideBy count 
                            |> Seq.roundTo !precision
                            |> Seq.toArray
                        subreadsByZmw 
                        |> Map.toSeq 
                        |> Seq.map fst 
                        |> Seq.toList 
                        |> List.sortWith Comparer.zmwIds
                        |> List.map (fun zmw -> (zmw, zmwWeight zmw))

                    // Write the calculated ZMW mapping scores to file
                    this.log Info "Writing per-zmw mapping scores to CSV..."
                    writeMappingCsv (prefix + "_zmws.csv") "ZmwId" zmwWeights
                // End of writeMappings block

                writePerBaseCsvResults "amplicon_analysis" results
                writePerReadCsvResults "amplicon_analysis_summary" results
                writeMappings "amplicon_analysis" results
            // End of writeConsensusResults enclosure

            writeResults results

        int ProcessExitCode.Success
