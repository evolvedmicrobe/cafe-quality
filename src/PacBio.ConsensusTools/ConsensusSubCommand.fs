namespace PacBio.ConsensusTools

open System
open PacBio.FSharp.Utils.ActivePatterns
open PacBio.FSharp.Utils.SubCommands
open PacBio.Consensus
open PacBio.Data
open ManyConsole
open ConsensusCore

(*
    A ConsensusSubCommand is a specialization of SubCommand that sets up 
    things needed for using the Quiver consensus algorithm.  Specifically:

    - Does things needed to initialize ConsensusCore
    - Sets up the chemistry parameters and allows override from the command line

 *)
[<AbstractClass>]
type ConsensusSubCommand(name, desc, ?hidden : bool) as this =
    inherit SubCommand(name, desc, defaultArg hidden false)

    override this.OverrideAfterHandlingArgumentsBeforeRun args =
        if this.Verbosity > 1 then
            ConsensusCore.Logging.EnableDiagnosticLogging()
        base.OverrideAfterHandlingArgumentsBeforeRun args

