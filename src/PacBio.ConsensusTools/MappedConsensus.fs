namespace PacBio.ConsensusTools

open System
open System.IO
open System.Linq
open System.Text
open System.Collections.Generic
open FSharpx.Collections

open PacBio.Consensus
open PacBio.Data
open PacBio.FSharp.Utils
open PacBio.FSharp.Utils.ActivePatterns
open PacBio.FSharp.Utils.Logging
open PacBio.FSharp.Utils.SubCommands
open PacBio.Utils


// --- Everything in these modules is going to be staged into PacBio.Consensus.FSharp
//     (with suitable renaming).

type CursorAlignedRead = IAlignedRead<IFullFeatures>
type ManifestAlignedRead = AlignedRead<QvFeatures>

type CursorRead = IRead<IFullFeatures>
type ManifestRead = Read<QvFeatures>

module CCHelper =

    let extractFeatures (read : #IRead<#IQvFeatures>) =
        let f = read.Features
        new ConsensusCore.QvSequenceFeatures(
            f.Basecalls(),
            f.InsertionQV(),
            f.SubstitutionQV(),
            f.DeletionQV(),
            f.DeletionTag(),
            f.MergeQV())

    let extractMappedRead origin (read : #IAlignedRead<#IQvFeatures>) =
        // TODO: assert aln is clipped
        let features = extractFeatures read
        let tStart = read.ReferenceInterval.Begin
        let tEnd   = read.ReferenceInterval.End
        let tStrand = enum<ConsensusCore.StrandEnum> ((int)read.Strand)
        let shortReadName = string <| read.LocationInFile
        let chem = match read.SequencingChemistry with
                   | null -> "unknown"
                   | x    -> x
        let ccRead = new ConsensusCore.Read(features, shortReadName, chem)
        new ConsensusCore.MappedRead(ccRead, tStrand, tStart - origin, tEnd - origin)


module Variants =

    type Variant = Variant

open Variants

module Consensus =

    open CCHelper

    let zerosLike (a : string) =
        Array.zeroCreate a.Length

    type Consensus =
        { Sequence   : string
          Confidence : int array
          Coverage   : int array }

        member this.Slice(beginIndex: int, len: int) =
          { Sequence   = this.Sequence.Substring(beginIndex, len);
            Confidence = this.Confidence.Slice(beginIndex, len);
            Coverage   = this.Coverage.Slice(beginIndex, len) }

        member this.Slice (window : Interval) =
            this.Slice(window.Begin, window.Length)

    type ErrorMessage = string

    type Call =
        | NoCall        of ErrorMessage
        | ConsensusCall of Consensus

    type NoCallStyle =
        | NoCallAsN
        | NoCallAsLowercaseReference

        static member Parse(s : string) =
            match s with
            | "N"                  -> NoCallAsN
            | "LowerCaseReference" -> NoCallAsLowercaseReference
            | _ -> invalidArg "s" "Unrecognized NoCallStyle requested"

    let noCallSequence style (refSeqInWin : string) =
        match style with
        | NoCallAsN -> String.replicate refSeqInWin.Length "N"
        | NoCallAsLowercaseReference -> refSeqInWin.ToLower()

    let intervals begin_ end_ windowLen  =
        seq { for b in begin_ .. windowLen .. (end_ - 1) do
                let e = min (b + windowLen) end_
                yield Interval(b, e) }

    let expandedWindow refLen expansion (window : Interval) =
        let b = max 0 (window.Begin - expansion)
        let e = min refLen (window.End + expansion)
        Interval(b, e)

    // equivalent to (sortByDescending fn |> Seq.truncate n), but does not require the
    // sortDescending to complete, because we know the optimum.
    let takeBest (n : int) (fn : 'a -> int) (optimum : int) (items : #seq<'a>)  : list<'a> =
        let itemsLL = LazyList.ofSeq items
        let fnCs = Func<'a, int>(fn)
        let rec go vals n (optimals : list<'a>) (subOptimals : list<'a>) =
            match (n, vals) with
            | (0, _) ->
                optimals
            | (_, LazyList.Nil) ->
                let rest = (subOptimals.OrderByDescending(fnCs)) |> Seq.truncate n |> Seq.toList
                List.append optimals rest
            | (_, LazyList.Cons(first, rest)) ->
                if (fn first) = optimum
                then
                    go rest (n-1) (first :: optimals) subOptimals
                else
                    go rest n optimals (first :: subOptimals)
        go itemsLL n List.empty List.empty


    let alnsForWindow (reader : IAlignmentReader) maxReads minMapQV refName win =
        // What would like to do is sort the (mapq-filtered) reads by template span within the
        // window and then take the top maxReads.  However if we are using a BAM file and the
        // coverage is very deep, it may be very slow to check the tStart, tEnd for all of them
        // (vs. cmp.h5 where it is very cheap).  Since in the common case we will find many reads that
        // completely span the query window immediately, we should have a means of quickly exiting,
        // without having to complete an entire sort and looking at all of the reads.  This
        // is what we get from the function `takeBest`, which implements a sort
        // that eagerly yields optimal values.
        let overlapWithWindow (aln : #IAlignedRead) = aln.ReferenceInterval.Intersection(win).Length

        reader.AlignmentsInInterval(refName, win)
        |> Seq.where (fun aln -> aln.MapQV > minMapQV)
        |> takeBest maxReads overlapWithWindow win.Length
        |> Seq.map AlignedRead<QvFeatures>.From
        |> Seq.map (fun a -> a.ReferenceClipped(win))
        |> Seq.toList


    // Given a set of reads clipped to the refWin, return the POA consensus sequence and
    // the Consensus.MappedReads orienting the reads with respect to the POA.
    // TODO: we should avoid using the reference---it should be easy to orient the reads using the
    // POA algorithm alone.
    let poaConsensusAndMappedReads (refWin : Interval) (refSeqInWin : string) (reads : list<#ManifestAlignedRead>) =

        let readIsFullLength (read : #ManifestAlignedRead) = read.ReferenceInterval.Spans(refWin)

        let readsForPoa =
            reads |> Seq.filter readIsFullLength
                  |> Seq.truncate 11
                  |> Seq.toArray

        let poa = if readsForPoa.Length < 3
                    then failwith "Too few reads for POA"
                    else
                        let sv = new ConsensusCore.StringVector(readsForPoa |> Array.map (fun r -> r.BasecallsForwardStrand()))
                        ConsensusCore.PoaConsensus.FindConsensus(sv)

        let poaCss = poa.Sequence()

        let poaToRefAln = ConsensusCore.ConsensusCore.Align(refSeqInWin, poaCss)
        let queryPositions = ConsensusCore.ConsensusCore.TargetToQueryPositions(poaToRefAln)

        // transform MappedRead into new coordinate space
        let lifted (mr : ConsensusCore.MappedRead) =
            let newStart = queryPositions.[mr.TemplateStart]
            let newEnd   = queryPositions.[mr.TemplateEnd]
            let copy = new ConsensusCore.MappedRead(mr)
            copy.TemplateStart <- newStart
            copy.TemplateEnd <- newEnd
            copy

        let poaMappedReads =  reads |> Seq.map (extractMappedRead refWin.Begin >> lifted) |> Seq.toList
        (poaCss, poaMappedReads)


    let consensusForMappedReads configTbl (refWin : Interval) (reads : list<ConsensusCore.MappedRead>) cssGuess =
        let mms = new ConsensusCore.SparseSseQvMultiReadMutationScorer(configTbl, cssGuess)
        for read in reads do
            mms.AddRead read |> ignore

        let converged = ConsensusCore.ConsensusCore.RefineConsensus(mms) // gross.
        if not converged then failwith "Did not converge to MLE" else ()

        ConsensusCore.ConsensusCore.RefineDinucleotideRepeats(mms)

        let tpl = mms.Template()
        let qvs = ConsensusCore.ConsensusCore.ConsensusQVs(mms).ToArray()
        let cov = zerosLike tpl

        { Sequence=tpl; Confidence=qvs; Coverage=cov; }

    let consensusForReads configTbl (refContig : ReferenceRecord) (refWin : Interval) (reads : list<ManifestAlignedRead>) =
        try
            let refSeqInWin = refContig.Sequence.Substring(refWin)
            let poaCss, mappedReads = poaConsensusAndMappedReads refWin refSeqInWin reads
            let css = consensusForMappedReads configTbl refWin mappedReads poaCss
            ConsensusCall css
        with
        | Failure(msg) -> NoCall msg

    let consensusForReads2 configTbl (refContig : ReferenceRecord) (refWin : Interval) (refDomain : Interval) (reads : list<ManifestAlignedRead>) =
        try
            let refSeqInWin = refContig.Sequence.Substring(refWin)
            let poaCss, mappedReads = poaConsensusAndMappedReads refWin refSeqInWin reads
            let css = consensusForMappedReads configTbl refWin mappedReads poaCss

            // Align to reference to identify correct clipping to allow windows to stich well.
            let ga = ConsensusCore.ConsensusCore.Align(refSeqInWin, css.Sequence)
            let targetPositions = ConsensusCore.ConsensusCore.TargetToQueryPositions(ga) |> Seq.toArray
            let cssStart = targetPositions.[refDomain.Begin - refWin.Begin]
            let cssLen   = targetPositions.[refDomain.End - refWin.Begin] - cssStart
            let css' = css.Slice(cssStart, cssLen)
            ConsensusCall css'
        with
        | Failure(msg) -> NoCall msg


open Consensus

// ----

module MappedConsensus =

    // TODO: hoist Window type up into PacBio.Data
    type Window(refName: string, interval: Interval) =
        member this.RefName = refName
        member this.Interval = interval
        override this.ToString() =  sprintf "%s:%d-%d" this.RefName this.Interval.Begin this.Interval.End
        static member Parse(winString) =
            match winString with
            | Regex "(.*):(\d+)-(\d+)" [refName; Integer winStart; Integer winEnd]
                -> Window (refName, new Interval(winStart, winEnd))
            | _ -> invalidArg "winString" "Invalid window string"

    type Focus =
        | FocusWindow of Window
        | FocusContig of string

        static member Parse(focusString) =
            match focusString with
            | Regex "(.*):(\d+)-(\d+)" _
                -> FocusWindow (Window.Parse(focusString))
            | Regex "([^:]*)" [refName]
                -> FocusContig refName
            | _ -> invalidArg "focusString" "Invalid focus string"

        static member ParseMultiple(focusString : string) =
            let subStrings = focusString.Split([|','|], StringSplitOptions.RemoveEmptyEntries) |> Seq.toList
            List.map Focus.Parse subStrings

        static member ParseMultipleFromFile(filename : string) =
            let lines = File.ReadAllLines(filename)
            let joinedLines = String.Join(",", lines)
            Focus.ParseMultiple(joinedLines)

        override this.ToString() =
            match this with
            | FocusWindow window  -> window.ToString()
            | FocusContig refName -> refName


    type ContigResult = { Consensus     : Consensus;
                          Variants      : Variant list }

    type Chunk =
        { Window  : Interval;  // Windows may be overlapping.
          Domain  : Interval;  // Interior of the window, non-overlapping; this is where our css and variants will pertain.
          Reads   : list<ManifestAlignedRead>  }

    type ChunkResult =
        { Domain        : Interval;
          Consensus     : Consensus;
          Variants      : Variant list }

    type MappedConsensus() as this =
        inherit ConsensusSubCommand ("MappedConsensus", "Consensus and variants from reads mapped to a scaffold or reference")

        let alnFilePath = ref ""
        let referencePath = ref ""
        let outputPath = ref "."
        let outputSuffix = ref None : (string Option) ref
        let refChunkSize = ref 500
        let refChunkOverlap = ref 5
        let noCallStyle = ref NoCallAsLowercaseReference
        let coverageDepth = ref 100
        let minMapQV = ref 10
        let foci = ref [] : (Focus list) ref

        do
            this.reqA "r|reference=" referencePath "Reference FASTA file."
            this.optF "f|foci="
                (fun f -> foci := Focus.ParseMultiple f)
                "Focal windows of the reference to consider (one or more comma-delimited strings like refName[:start-end])"
            this.optF "F|fociFile="
                (fun f -> foci := Focus.ParseMultipleFromFile f)
                "Process (newline delimited) focal windows specifed in the given file"
            this.optA "m|minMapQV="  minMapQV "Minimum mapping QV"
            this.optA "o|output=" outputPath "Directory to write results files. Will be created if it doesn't exist."
            this.optF "s|outputFilenameSuffix="
                (fun arg -> outputSuffix := Some arg)
                "Suffix for output filenames; will output files like consensus.$SUFFIX.fasta"
            this.optF "noCallStyle="
                (fun ncs -> noCallStyle := NoCallStyle.Parse(ncs))
                "Output style for nocalled bases/windows: either 'N' or 'LowerCaseReference' (default)"

            this.AllowsAnyAdditionalArguments () |> ignore

        let report (refContig : ReferenceRecord) (fastaSink : FastaWriter) gffSink  (r : ContigResult) =
            let cssContigName = refContig.Name + "|quiver"
            fastaSink.WriteRecord(cssContigName, r.Consensus.Sequence)

        let flattenNoCalls noCallStyle (refSeqInWin : string) (c : Call) =
            match c with
            | NoCall _   ->  { Sequence   = noCallSequence noCallStyle refSeqInWin;
                               Confidence = zerosLike refSeqInWin;
                               Coverage   = zerosLike refSeqInWin }
            | ConsensusCall cc -> cc

        let processChunk configTbl (refContig : ReferenceRecord) (chunk : Chunk) : ChunkResult =
            let css = consensusForReads2 configTbl refContig chunk.Window chunk.Domain chunk.Reads
            match css with
                | ConsensusCall _ -> this.logf Info "Called consensus for %O" chunk.Domain
                | NoCall reason   -> this.logf Info "Nocall for %O (reason: %s)" chunk.Domain reason
            let refSeqInDomain = refContig.Sequence.Substring(chunk.Domain)
            { Domain    = chunk.Domain;
              Consensus = flattenNoCalls !noCallStyle refSeqInDomain css;
              Variants  = List.empty }

        let processContig configTbl (contig : ReferenceRecord) (focus : Interval option) reader : ContigResult =
            let wins = match focus with
                        | Some focusInterval -> intervals focusInterval.Begin focusInterval.End !refChunkSize
                        | None -> intervals 0 contig.Length !refChunkSize

            let chunks = seq { for domain in wins ->
                                let exWin = expandedWindow contig.Length !refChunkOverlap domain
                                { Window = exWin
                                  Domain = domain;
                                  Reads = alnsForWindow reader !coverageDepth !minMapQV contig.Name exWin; } }
            let resultChunks = Seq.parMap' (processChunk configTbl contig) chunks

            let cssSequence = new StringBuilder()
            let cssConfidence = new List<int>()
            let cssCoverage = new List<int>()
            let contigVariants = new List<Variant>()

            for chunk in resultChunks do
                let css = chunk.Consensus
                cssSequence.Append(css.Sequence) |> ignore
                cssConfidence.AddRange(css.Confidence) |> ignore
                cssCoverage.AddRange(css.Coverage) |> ignore
                contigVariants.AddRange(chunk.Variants)

            let contigCss = { Sequence   = cssSequence.ToString();
                              Confidence = cssConfidence.ToArray();
                              Coverage   = cssCoverage.ToArray(); }
            { Consensus = contigCss;
              Variants  = Seq.toList contigVariants }



        override this.Run(remainingArgs) =
            alnFilePath :=
                if remainingArgs.Length >= 1
                    then Seq.head remainingArgs
                    else invalidOpt "Alignment file must be specified"

            use reader = AlignmentReader.Open(!alnFilePath)
            let scorerConfig = this.LoadQuiverConfig "QuiverParameters.ini" (set reader.PulseFeatures) (set reader.SequencingChemistries)

            this.PrepareOutputDirectory !outputPath
            let consensusFastaFname =
                match !outputSuffix with
                | None ->        "consensus.fasta"
                | Some suffix -> "consensus." + suffix + ".fasta"
            use consensusFasta = new FastaWriter(Path.Combine(!outputPath, consensusFastaFname))
            reader.ReferenceTable.LoadSequencesFromFasta(!referencePath)
            let variantsGff = ()

            let configTbl' = scorerConfig.Parameters
            // HACK: LAA and CCS use ridiculously large banding.  Hack in narrower banding here until
            // we can switch all code over to narrower bands.
            let configTbl = new ConsensusCore.QuiverConfigTable()
            for key in configTbl'.Keys() do
                let config = configTbl'.At(key)
                if (not scorerConfig.HasChemistryOverride) && (key = "*")
                    then
                        configTbl.Insert(key, config) |> ignore
                    else
                        let adjustedConfig = new ConsensusCore.QuiverConfig(config)
                        adjustedConfig.Banding.ScoreDiff <- 6.0f
                        configTbl.Insert(key, adjustedConfig) |> ignore

            let processFocus focus =
                try
                    match focus with
                    | FocusWindow window ->
                        let focusContig = reader.ReferenceTable.Lookup(window.RefName)
                        this.logf Info "Processing window %O" window
                        processContig configTbl focusContig (Some window.Interval) reader |> report focusContig consensusFasta variantsGff
                        this.logf Info "Finished window %O" window
                    | FocusContig refName ->
                        let focusContig = reader.ReferenceTable.Lookup(refName)
                        if reader.HasAlignmentsForReference(refName)
                        then
                            this.logf Info "Processing contig %O" refName
                            processContig configTbl focusContig None reader |> report focusContig consensusFasta variantsGff
                            this.logf Info "Finished contig %O" refName
                        else
                            ()
                with :? KeyNotFoundException ->
                    this.logf Info "Skipping focus '%O' with unrecognized contig" focus


            match !foci with
                | [] ->
                    let foci' = seq { for contig in reader.ReferenceTable -> FocusContig contig.Name }
                    foci' |> Seq.iter processFocus
                | _ ->
                    !foci |> Seq.iter processFocus

            int ProcessExitCode.Success


    open NUnit.Framework

    [<TestFixture>]
    type TestTakeBest() =

        [<Test>]
        member this.Test1() =
            let vals = Seq.concat [ seq { 0..99 };
                                    Seq.init 10 (fun _ -> 100);
                                    seq { for i in 1..10 -> i / 0 } ]
            let ans = takeBest 10 id 100 vals
            Assert.AreEqual(List.init 10 (fun _ -> 100), ans)
            ()
