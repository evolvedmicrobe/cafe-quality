
namespace PacBio.CCSTrain

open System.IO

open PacBio.Consensus
open PacBio.FSharp.Utils.Logging
open PacBio.FSharp.Utils.SubCommands
open PacBio.Utils

type CCSTrainCompare() as this =
    inherit SubCommand ("CCSTrainCompare", "Train CCS", hidden = true)

    let reference = ref ""
    let alignments = ref ""
    let fofn = ref ""
    let maxTraces = ref 30
    let param1 = ref ""
    let param2 = ref ""

    do
        this.optA "r|ref=" reference "Reference FASTA file."
        this.optA "l|alignments=" alignments "Alignment cmp.h5 file."
        this.optA "f|fofn=" fofn "File of filenames for bas.h5."
        this.optA "t|maxTraces=" maxTraces "Maximum number of traces to process."
        this.optA "param1=" param1 "First Parameter Set To Use"
        this.optA "param2=" param2 "Second Parameter Set To Use"
    
    override this.Run(args) =
        let tester = new TrainCompare ()
        tester.Run (!alignments, !fofn, !reference, !maxTraces, !param1, !param2) |> ignore
        0