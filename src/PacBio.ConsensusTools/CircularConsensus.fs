namespace PacBio.ConsensusTools

open System.IO

open PacBio.Consensus
open PacBio.IO
open PacBio.FSharp.Utils
open PacBio.FSharp.Utils.ActivePatterns
open PacBio.FSharp.Utils.Logging
open PacBio.FSharp.Utils.SubCommands
open PacBio.Utils

type CircularConsensus() =

    let fofn = ref ""
    let zmwStart = ref 0
    let zmwCount = ref 0

 
    member this.Run (args : string[]) =
        // CircularConsensus -n 8 -o test -fofn=/Users/nigel/CCS_P6_C4/input.fofn  -csv
        fofn := args.[0]

       
        let processMoviePart basFile =
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

        let basFiles = Fofn.Filenames !fofn
           
        for basFile in basFiles do
            let path = Path.GetFileName basFile
            let parts = path.Split '.'
            let stem =
                if System.Int32.TryParse (parts.[1], ref 0) then
                    String.concat "." (Seq.take 2 parts)
                else
                    parts.[0]
            
            
            processMoviePart basFile

        int ProcessExitCode.Success
