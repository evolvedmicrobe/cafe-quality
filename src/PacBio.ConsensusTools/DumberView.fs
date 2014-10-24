
namespace PacBio.ConsensusTools

open PacBio.FSharp.Utils.SubCommands

type DumberView() =
    inherit SubCommand("DumberView", "", hidden = true)
    override this.Run(args) = 0
