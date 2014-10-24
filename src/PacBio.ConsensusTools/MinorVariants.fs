
namespace PacBio.ConsensusTools

open System.Linq
open System.IO

open PacBio.Consensus
open PacBio.IO
open PacBio.IO.Fasta
open PacBio.FSharp.Utils.ActivePatterns
open PacBio.FSharp.Utils.Logging
open PacBio.FSharp.Utils.SubCommands
open PacBio.Utils

type MinorVariants() as this = 
    inherit SubCommand ("MinorVariants", "Detect and characterize minor populations in long aligned reads")

    let cmpH5Path = ref ""
    let referencePath = ref ""
    let outputPath = ref "."
    let maxMolecules = ref 0
    let refChunkSize = ref 2000
    let minLength = ref 50
    let minAccuracy = ref 0.9
    let minLocalAccuracy = ref 0.9
    let minErrorRate = ref 0.005
    let minCoverage = ref 100
    let windowRadius = ref 10
    let maxPValue = ref 0.05
    let confidenceInterval = ref 0.95
    let doAminoVariants = ref false
    let doCountAmbiguous = ref false

    do
        this.reqA "r|reference=" referencePath "Reference FASTA file."
        this.reqA "cmp=" cmpH5Path "Input cmp.h5 file."
        this.optA "o|output=" outputPath "Directory to write results files. Will be created if it doesn't exist."
        this.optA "maxMolecules=" maxMolecules "Maximum number of molecules to include in analysis. Pass 0 for use everything available."
        this.optA "referenceChunkSize=" refChunkSize "Maximum size of a reference chunk to process at a time."
        this.optA "minLength=" minLength "Minimum subread length to use."
        this.optA "minAccuracy=" minAccuracy "Minimum subread alignment accuracy."
        this.optA "minLocalAccuracy=" minLocalAccuracy "Minimum local re-alignment accuracy."
        this.optA "minErrorRate=" minErrorRate "The floor of the minimum estimated error rate."
        this.optA "minCoverage=" minCoverage "Minimum coverage required to score a variant."
        this.optA "windowRadius=" windowRadius "Radius of window to use when doing local re-alignment."
        this.optA "maxPValue=" maxPValue "Maximum P-value to report a variant."
        this.optA "confidenceInterval=" confidenceInterval "Confidence interval about variant frequency."
        this.optF "aminoVariants" (fun _ -> doAminoVariants := true) "Translate reference and call variants in amino acid space. Default = false"
        this.optF "countAmbiguous" (fun _ -> doCountAmbiguous := true) "When sequences at a site map ambiguously, count all possibilities rather than none. Default = false"

        this.AllowsAnyAdditionalArguments () |> ignore

    override this.Run(args) =

        let logger = Logger "MinorVariants"
        let log lvl msg = logger.Log lvl msg
        let logf lvl fmt = logger.LogFormat lvl fmt

        log Info (System.Environment.GetCommandLineArgs () |> String.concat " ")

        if File.Exists !outputPath then
            invalidOptf "Ouput path: '%s' already exists as a file" !outputPath

        if Directory.Exists !outputPath then
            logf Debug "Writing to existing results folder: '%s'" !outputPath
        else
            Directory.CreateDirectory !outputPath |> ignore

        if !confidenceInterval < 0.0 || !confidenceInterval > 1.0 then
            invalidOptf "Confidence interval must be in [0, 1], got '%g'" !confidenceInterval
        else
            let cmpH5 =
                let reader =
                    match (!cmpH5Path, args) with
                    | Named(path) ->
                        if args.Length > 0 then
                            invalidOptf "Unknown arguments: '%s'" <| String.concat " " args
                        else
                            OldCmpH5Reader.CreateReader path
                    | Unnamed(paths) -> OldCmpH5Reader.CreateReader paths.[0]
                    | NonExistent(paths) ->
                        invalidOptf "Input file(s) '%s' do(es) not exist -- specify .bas or .bax files on the command line or use the --fofn argument" (String.concat "', '" paths)
                    | _ ->
                        invalidOpt "Input files unspecified -- specify .bas or .bax files on the command line or use the --fofn argument"
                if reader.IsSorted then
                    reader
                else
                    invalidOpt "Input cmp.h5 is not sorted"

            let references =
                match !referencePath with
                | File(path) ->
                    let reader = new SimpleFASTAReader (!referencePath)
                    seq { for contig in reader do yield (contig.Header, contig.GetSequence ()) } |> Map.ofSeq
                | _ ->
                    invalidOptf "Reference file '%s' is empty or does not exist" !referencePath

            let config = VariantDetectionConfig(
                             CountAmbiguous = !doCountAmbiguous,
                             WindowRadius = !windowRadius,
                             MinLocalAccuracy = !minLocalAccuracy,
                             ConfidenceInterval = !confidenceInterval,
                             ConvergenceCriterion = 1E-8,
                             MaxIterations = 1000,
                             MaxPValue = !maxPValue,
                             MinCoverage = !minCoverage,
                             BinomialTestThreshold = 0.01,
                             MinErrorRate = !minErrorRate)

            let refPath = Path.GetFileName !referencePath
            use csvWriter = new CsvVariantWriter (Path.Combine (!outputPath, "minor_variants.csv"), refPath)
            use vcfWriter = new VcfWriter (Path.Combine (!outputPath, "minor_variants.vcf"), refPath)

            for contig in references do
                let fastaId = contig.Key
                let sequence = contig.Value

                let refInfo =
                    let ri = cmpH5.References.First (fun r -> r.FullName = fastaId)
                    match ri with
                    | null -> 
                        invalidOptf "Couldn't find contig '%s' in cmp.h5 -- are you using the correct reference FASTA?" fastaId
                    | _ -> ri
                
                let doTranslate =
                    try
                        do Translation.Translate sequence |> ignore
                        !doAminoVariants
                    with
                    | _ -> false

                let csvContig = csvWriter.AddContig (fastaId, sequence, doTranslate)
                let vcfContig = vcfWriter.AddContig (fastaId, sequence, doTranslate)

                let isValid = fun (aln : IAlnSummary) -> aln.ReadLength >= !minLength ||
                                                         aln.Accuracy >= !minAccuracy

                for start in 0 .. !refChunkSize .. (refInfo.Length - 1) do
                    let end_ = min refInfo.Length (start + !refChunkSize)
                    logf Info "Loading reads in range [%i, %i]" (start + 1) end_
                    let alns =
                        let all = [| for aln in cmpH5.GetReadsInRange(refInfo.ID, start, end_) do
                                     if isValid aln then yield aln |]
                        if !maxMolecules > 0 then
                            all.ReservoirSample !maxMolecules
                        else
                            all
                    if alns.Length = 0 then
                        logf Info "No valid alignments found in range [%i, %i]" (start + 1) end_
                    else
                        logf Info "Filtered %i alignments, detecting variants" alns.Length
                        let chunk = sequence.Substring (start, end_ - start)
                        let ttpl = new TrialTemplate (chunk)
                        let detect = new VariantDetection (config)
                        let variants = detect.LearnApply (alns, ttpl, start) |> Seq.toList

                        csvContig.AddVariants variants
                        vcfContig.AddVariants variants

            int ProcessExitCode.Success
