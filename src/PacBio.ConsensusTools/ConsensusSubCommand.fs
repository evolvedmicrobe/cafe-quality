namespace PacBio.ConsensusTools

open System
open PacBio.FSharp.Utils.ActivePatterns
open PacBio.FSharp.Utils.SubCommands
open PacBio.Consensus
open PacBio.Data

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

    let chemistry = ref ""
    let parameters = ref ""

    do
        this.optA "p|parameters=" parameters "Algorithm parameters path (defaults to built-in parameters)"
        this.optA "c|chemistry=" chemistry "Use a specified chemistry model. Default = autodetection"

    override this.OverrideAfterHandlingArgumentsBeforeRun args =
        if this.Verbosity > 1 then
            ConsensusCore.Logging.EnableDiagnosticLogging()
        base.OverrideAfterHandlingArgumentsBeforeRun args

    member this.ChemistryModelOverride with get () = !chemistry

    member this.LoadQuiverConfig (chemistryParameterFilebase : string) (pulseFeaturesInFile : string Set) (chemistriesInFile : string Set) =
        let paramsFile = ParameterLoading.SelectParameterFile (!parameters, chemistryParameterFilebase)
        match !chemistry with
        | ChemistryMapping(_)         -> invalidOpt "Passing chemistry information by chemistry_mapping.xml is no longer supported."
        | ChemistryModel(chem, model) -> ParameterLoading.LoadParametersFromFile (paramsFile, chem, model)
        | NonEmpty(chem)              -> ParameterLoading.LoadParametersFromFile (paramsFile, chem)
        | _                           ->
            // Chemistry info must be available in the data file(s)
            // FIXME: check to make sure that all the chemistries asked for are loaded from the parameters
            if chemistriesInFile.IsEmpty
            then raise (ChemistryLookupException "Chemistry information absent from data file(s)")
            else ParameterLoading.LoadParametersFromFile paramsFile