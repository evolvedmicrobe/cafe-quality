
namespace PacBio.CCSTrain

open System.IO

open PacBio.Consensus
open PacBio.FSharp.Utils.Logging
open PacBio.FSharp.Utils.SubCommands
open PacBio.Utils

type CCSTrain() as this =
    inherit SubCommand ("CCSTrain", "Train CCS", hidden = true)

    let reference = ref ""
    let alignments = ref ""
    let fofn = ref ""
    let outputFile = ref ""
    let maxIterations = ref 1000
    let maxTraces = ref 30
    let snrCut = ref SnrCut.PassAll

    do
        this.optA "r|ref=" reference "Reference FASTA file."
        this.optA "l|alignments=" alignments "Alignment cmp.h5 file."
        this.optA "f|fofn=" fofn "File of filenames for bas.h5."
        this.optA "o|output=" outputFile "Root output filename."
        this.optA "i|maxIterations=" maxIterations "Maximum number of iterations. Pass -1 for no maximum."
        this.optA "t|maxTraces=" maxTraces "Maximum number of traces to process."
        this.optF "snrCut=" (fun s -> snrCut := SnrCut.Parse(s)) "Only process ZMWs within a given window of SNR space"

    override this.Run(args) =
        let trainer = new ConsensusTrain ()

        match trainer.Run (!alignments, !fofn, !reference, !outputFile, !maxTraces, !maxIterations, !snrCut) with
        | 0 -> int ProcessExitCode.Success
        | _ -> int ProcessExitCode.UnspecifiedError
    