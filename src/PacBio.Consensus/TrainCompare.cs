using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Diagnostics;
using System.Linq;
using ConsensusCore;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics.LinearAlgebra;

using PacBio.Utils;

using Vector=MathNet.Numerics.LinearAlgebra.Double.DenseVector;
namespace PacBio.Consensus
{
    /// <summary>
    /// Class that compares the objective function and accuracy for two different sets of parameters on the same training function.
    /// </summary>
    public class TrainCompare
    {

            private static PacBioLogger logger = PacBioLogger.GetLogger("TrainCompare");

            private static void Log(string msg)
            {
                logger.Log(LogLevel.INFO, msg);
            }

            private static void Log(LogLevel level, string msg)
            {
                logger.Log(level, msg);
            }

            private static void Log(LogLevel level, string msg, params object[] args)
            {
                logger.Log(level, String.Format(msg, args));
            }

            private static void Log(string msg, params object[] args)
            {
                logger.Log(LogLevel.WARN, String.Format(msg, args));
            }

        public void TestParameters(List<CCSExample> data, string params1, string params2 )
        {           
                
            foreach (string ps in new string[] {params1, params2}) {

                using (var scConfig = ParameterLoading.DefaultCCS)
                using (var qvConfig = scConfig.Parameters.At (ps))
                using (var start = qvConfig.QvParams) {
                    var algo = RecursionAlgo.Viterbi;

                    // Train the model
                    float trainError;
                    float testError;
                    Console.WriteLine ("Results for: " + ps);
                    var result = LikelihoodObjective (qvConfig.QvParams, data, algo);
                    Console.WriteLine ("Score is: " + result);
                    // Assess the accuracy on test set
                    Log ("Measuring Accuracy:");
                    AccuracyObjective (qvConfig.QvParams, data, algo);
                }
            }

        }



        /// <summary>
        /// Subcommands must implement this method to carry out the subcommand.
        /// </summary>
        public void Run(string cmpH5File, string fofn, string reference, int totalTraces, string param1, string param2)
        {                

            Log (LogLevel.WARN, "Starting to optimize across partitions");

            // Load reference
            var rFasta = new PacBio.IO.Fasta.SimpleFASTAReader(reference);
            var refDict = rFasta.ToDictionary(e => e.Header, e => e.GetSequence());

            var totalReferences = refDict.Count;
            Log (LogLevel.WARN, "Found " + totalReferences + " references in file");
            var tracesPerReference = (int) Math.Ceiling (totalTraces / (double)totalReferences); 


            List<CCSExample> data = new List<CCSExample> (totalTraces);
            foreach (string refer in refDict.Keys) {
                var traceSet = new TraceSet (cmpH5File, fofn);

                var cur = CCSExample.GetSimpleExamples (traceSet, refDict, tracesPerReference, refer);
                data.AddRange (cur);
            }
            TestParameters (data, param1, param2);

        }




            public float LikelihoodObjective(QvModelParams spec, List<CCSExample> data, RecursionAlgo algo)
            {
                using (var scConfig = new ScorerConfig { Parameters = new QuiverConfigTable(), Algorithm = algo, HasChemistryOverride = true })
                using (var qvConfig = new QuiverConfig(spec, (int) Move.ALL_MOVES, new BandingOptions(4, 18), -12.5f))
                {
                    scConfig.Parameters.Insert("*", qvConfig);

                    // FIXME (lhepler) -- pipe through the number of jobs
                    var diffsAndOverall = data.ParSelect(
                        e =>
                        {
                            // Evaluate the likelihood of the correct answer under a logit model
                            var trueTpl = e.CorrectTrialTemplate;
                            var bases = e.Trace.ZmwBases;

                            using (var scorer = new MultiReadMutationScorer(e.Regions, bases, e.CorrectTrialTemplate, scConfig))
                            {
                                var overallScore = scorer.GetBaselineScores().Sum();

                                var muts = GenerateMutations.GenerateUniqueMutations(trueTpl);

                                var ll = muts.Select(
                                    m =>
                                    {
                                        var score = scorer.ScoreMutation(m);
                                        return ( 1 / (1 + Math.Exp(-score)));
                                    }).Sum();

                                return new Tuple<float, float>((float) ll, overallScore);
                            }
                        }).ToArray();

                    var r = diffsAndOverall.Map(v => v.Item1);
                    var overall = diffsAndOverall.Map(v => v.Item2);

                    //Log("Current model:");
                    //Log(spec.ToString());

                    var maxErr = r.Max();
                    var minErr = r.Min();
                    var mean = r.Average();

                    Log(LogLevel.WARN, "Mean LL: {0:0.0000}, Max: {1:0.0000}, Min: {2:0.0000}, Mean Aln Score: {3:0.00}", mean, maxErr, minErr, overall.Sum());

                    Log(LogLevel.INFO, "Parameters:\n{0}", 
                        Vector.OfEnumerable(ConsensusCoreWrap.QvModelParamsToArray(spec).Select(v => (double) v)).ToString());

                    return mean;
                }   
            }

            private CCSStream ccsAlgo = CCSStream.TrainingConfig;
            public float AccuracyObjective(QvModelParams spec, List<CCSExample> data, RecursionAlgo algo)
            {
                using (var scConfig = new ScorerConfig { Parameters = new QuiverConfigTable(), Algorithm = algo, HasChemistryOverride = true })
                using (var qvConfig = new QuiverConfig(spec, (int) Move.ALL_MOVES, new BandingOptions(4, 18), -12.5f))
                {
                    scConfig.Parameters.Insert("*", qvConfig);

                    ccsAlgo.ScConfig = scConfig;

                    var rStart = data.ParSelect(
                        e =>
                        {
                            var newConsensus = ccsAlgo.Map(e.Trace.ZmwBases).Item2;

                            var al = CCSExample.LocalAlign(
                                e.Reference,
                                newConsensus.Sequence);

                            if (al.MiscallRate < 0)
                            {
                                Log("Got bad result: {0} - NBases: {1}", al.Alignment, newConsensus.NumBases);
                                if (Debugger.IsAttached)
                                {
                                    Debugger.Break();
                                }
                            }

                            return al;

                        }).ToArray();

                    // the scConfig will be disposed at the end of this block
                    ccsAlgo.ScConfig = null;

                    Log("Current model:");
                    Log(spec.ToString());

                    var nEmpty = rStart.Where(a => a.TemplateLength < 1).Count();

                    // Cut out bad CCS reads -- they can impact the accuracy we see
                    var r = rStart.Where(a => a.TemplateLength > 1).ToArray();

                    if (nEmpty > 0)
                        Log("Got {0} empty alignments.", nEmpty);

                    var ir = r.Map(a => a.InsertionRate).Average();
                    var dr = r.Map(a => a.DeletionRate).Average();
                    var mr = r.Map(a => a.MiscallRate).Average();
                    var er = r.Map(a => a.ErrorRate).Average();

                    var maxErr = r.Max(a => a.ErrorRate);
                    var minErr = r.Min(a => a.ErrorRate);

                    Log("Insertion Rate: {0:0.0000}, Deletion Rate: {1:0.0000}, Miscall Rate: {2:0.0000}, Error Rate: {3:0.0000}, Max: {4:0.0000}, Min: {5:0.0000}",
                        ir, dr, mr, er, maxErr, minErr);

                    // Use a robust measure of the accuracy, so that we don't train for outliers
                    //var final = r.Average(v => 1f - (float)v.Accuracy);

                    // We will seek the 'minimax' error rate, to try and balance the error modes
                    var final = (float)Math.Max(dr, Math.Max(ir, mr));
                    return final;
                }   
        }
    }
}

