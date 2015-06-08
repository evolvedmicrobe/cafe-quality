using System;
using System.Collections.Generic;
using PacBio.IO;
namespace PacBio.Consensus
{
    public static class RInterop
    {
        public static Dictionary<string, object> GetSubReadsForZMW(string basFile, double zmw) 
        {
            var basReader = BaseReader.CreateSource (basFile);
            var chemistries = basReader.SequencingChemistry;

            var config = new ConsensusConfig() {
                MinPredictedAccuracy = 0.0f,
                MinFullPasses = 0,
                MinLength = 0,
                MaxLength = Int32.MaxValue,
                SnrCut = SnrCut.PassAll};
            
            var stream = new CCSStream (config);
            var bases = basReader.ByHoleNumber (Convert.ToInt32(zmw));
            return stream.ReturnAlignments (bases);

        }
        public static Dictionary<string, object> WTF()
        {
            var toR = new Dictionary<string, object> ();
            toR ["test"] = new double[4];
            toR ["jj"] = "ASDfadfa";
            return toR;
        }


//            let chemInfo = basReader.ChemistryBarcode
//            this.logf Info "Got SequencingChemistry '%s' from ba[sx].h5/metadata.xml" basReader.SequencingChemistry
//            h5Sink.AddChemistryInformation (chemInfo.BindingKit, chemInfo.SequencingKit, chemInfo.ChangelistID)
//
//            let sinkFunc (rawBases, consensusBasesResult) =
//            let consensusBases = snd consensusBasesResult
//            fastaSink.OnNext consensusBases
//            h5Sink.OnNext consensusBases
//            csvSink.OnNext (rawBases, consensusBases)
//            zmwResultReport.TallyResult((fst consensusBasesResult))
//            ()
//
//            try
//            toProcess |> Seq.parMap (fun x -> (x, stream.Map x)) |> Seq.iter sinkFunc
//            fastaSink.OnComplete()
//            h5Sink.OnComplete()
//            csvSink.OnComplete()
//            with
//            | ex ->
//            fastaSink.OnError(ex)
//            h5Sink.OnError(ex)
//            csvSink.OnError(ex)
//            reraise ()
//
//            let basFiles =
//                match (!fofn, args) with
//                | Named(path) ->
//                if args.Length > 0 then
//                    invalidOptf "Unknown arguments: '%s'" <| String.concat " " args
//                else
//                    Fofn.Filenames path
//                    | Unnamed(paths) -> paths
//                    | NonExistent(paths) ->
//                    invalidOptf "Input file(s) '%s' do(es) not exist -- specify .bas or .bax files on the command line or use the --fofn argument." (String.concat "', '" paths)
//                    | _ ->
//                    invalidOpt "Input files unspecified -- specify .bas or .bax files on the command line or use the --fofn argument."
//
//                    for basFile in basFiles do
//                        let path = Path.GetFileName basFile
//                        let parts = path.Split '.'
//                        let stem =
//                            if System.Int32.TryParse (parts.[1], ref 0) then
//                                String.concat "." (Seq.take 2 parts)
//                            else
//                                parts.[0]

        public static double RInteropCall ()
        {
            Console.WriteLine ("Call made");
            var l = new List<string>(2);
            l.Add("dump");
            l.Add("dump1");
            return 0.0;
        }
    }
}

