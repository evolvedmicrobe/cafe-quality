
namespace PacBio.ConsensusTools
open System
open PacBio.HDF
open PacBio.FSharp.Utils.Logging
open PacBio.FSharp.Utils.SubCommands
open PacBio.Utils
open PacBio.Data

module EntryPoint =


    [<EntryPoint>]
    let main args =

        let code = CircularConsensus()
        code.Run args |> ignore
        0