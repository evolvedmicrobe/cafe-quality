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
    // (more general implementation is coming in PacBio.Data)
    internal struct Range
    {
        public readonly float Begin;
        public readonly float End;

        public Range(float begin, float end) 
        {
            Begin = begin;
            End = end;
        }

        public static readonly Range Universe = new Range(float.MinValue, float.MaxValue);

        public bool Contains(float pt)
        {
            return (Begin <= pt && pt < End);
        }

        public static Range Parse(string s)
        {
            var elts = s.Split(new char[] {'-'});
            return new Range(float.Parse(elts[0]), float.Parse(elts[1]));
        }
    }
        
    public class SnrCut
    {
        // in TGAC order
        private Range[] Ranges;
        private SnrCut() {}

        public static readonly SnrCut PassAll = new SnrCut { Ranges =  4.Fill(Range.Universe) };

        private static Tuple<char, Range> MakeRange(string clause)
        {   
            char channel = clause[0];
            var range = Range.Parse(clause.Substring(1));
            return Tuple.Create(channel, range);
        }

        // format expected: A3-5,T6-8
        // omitted channels have no filtering applied
        public static SnrCut Parse(string s)
        {
            var clauses = s.Split(new char[] {','}, StringSplitOptions.RemoveEmptyEntries);
            var ranges = 4.Fill(Range.Universe);
            clauses.Select(MakeRange).ForEach(t =>
            {
                switch(t.Item1) {
                    case 'T': ranges[0] = t.Item2; break;
                    case 'G': ranges[1] = t.Item2; break;
                    case 'A': ranges[2] = t.Item2; break;
                    case 'C': ranges[3] = t.Item2; break;
                }
            });
            return new SnrCut { Ranges = ranges };
        }

        public bool Contains(float[] snrs)
        {
            return Ranges.Zip(snrs, (r, snr) => r.Contains(snr)).All();
        }
    }

    public class ConsensusTrain
    {		

        const float FIXED_MATCH_SCORE = 1.0F;

        private static PacBioLogger logger = PacBioLogger.GetLogger("ConsensusTrain");
		
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

        private IOptimizer optimizer;

        private ReadConfigurationAssigner rca = new ReadConfigurationAssigner ();

		public ConsensusTrain()
        {

        }

        /// <summary>
        /// Train a single snr/coverage group
        /// </summary>
        /// <returns>The run.</returns>
        public int InnerRun(List<CCSExample> train, List<CCSExample> test, int hopedForTraces,
            int snrGroup, int coverageGroup, string outFile, int maxIterations = 100 )
        {          
            try {
            // Check that enough data is available
            int warnThresh = hopedForTraces / 2;
            int failThresh = hopedForTraces / 10;

            if (train.Count < failThresh || test.Count < failThresh)
            {
                Log(LogLevel.ERROR, "Insufficient data for training: maxTraces={0}, train.Count={1}, test.Count={2}",
                    hopedForTraces, train.Count, test.Count);

                return 1;
            }

            if (train.Count < warnThresh || test.Count < warnThresh)
            {
                Log(LogLevel.WARN, "Input data yields less than 2x requested: maxTraces={0}, train.Count={1}, test.Count={2}",
                    hopedForTraces, train.Count, test.Count);
            }

            using (var scConfig = ParameterLoading.DefaultCCS)
            using (var qvConfig = scConfig.Parameters.At("P6-C4"))
            using (var start = qvConfig.QvParams)
            {
                var algo = RecursionAlgo.Prob;

                // Train the model
                float trainError;
                float testError;
                var trainedModel = OptimizeConsensus(train, test, start, algo, maxIterations, out trainError, out testError);

                // Assess the accuracy on test set
                Log("Measuring Accuracy:");
                AccuracyObjective(trainedModel, test, algo);

                // Save the model to a results file
                Log("Writing CCS MoveSpec to: {0}", outFile);

                // FIXME (lhepler) pass chemistry and model type appropriately
                var name = ReadConfigurationAssigner.GetNameForConfiguration (snrGroup, coverageGroup);
                ParameterLoading.SaveModelToFile(name, "AllQVsMergingByChannelModel", trainedModel, outFile);
            }
            return 0;
            }
            catch(Exception thrown) {
                Console.WriteLine ("Failed!!!\n\n" + thrown.Message + "\n\n" + thrown.StackTrace);
               
            }
            return 1;
        }



        /// <summary>
        /// Subcommands must implement this method to carry out the subcommand.
        /// </summary>
        public int Run(string cmpH5File, string fofn, string reference, string outFile, int totalTraces, 
            int maxIterations = 100, SnrCut snrCut = null)
        {
           
            if (String.IsNullOrEmpty(outFile))
            {
                Log(LogLevel.ERROR, "Specify an output file with '-o <filename>'");
                return 1;
            }
            if (System.IO.File.Exists (outFile)) {
                System.IO.File.Delete (outFile);
            }

            Log (LogLevel.WARN, "Starting to optimize across partitions");


            // Load reference
            var rFasta = new PacBio.IO.Fasta.SimpleFASTAReader(reference);
            var refDict = rFasta.ToDictionary(e => e.Header, e => e.GetSequence());

            var totalReferences = refDict.Count;
            Log (LogLevel.WARN, "Found " + totalReferences + " references in file");
            var tracesPerReference = (int) Math.Ceiling (totalTraces / (double)totalReferences); 


            // Load training data for each partition
            // This code goes through the entire file and assigns examples to different partitions
            var traceSet = new TraceSet(cmpH5File, fofn);
            var tds = CCSExample.GetExamples (traceSet, refDict, tracesPerReference * 2, rca);

            tds.PrintCounts ();

            // Now assign all the training sets to different files
            for (int i = 0; i < rca.NumberOfSNRGroups; i++) {
                for (int j = 0; j < rca.NumberOfCoverageGroups; j++) {
                    Log (LogLevel.WARN, "Optimizing " + i +" - "+ j);
                    // Run this data for this file and output it.
                    var traintest = tds.GetExamples (i, j, tracesPerReference);
                    var train = traintest.Item1;
                    var test = traintest.Item2;
                    var res = InnerRun (train,test, totalTraces, i, j,outFile, maxIterations);       
                    if (res != 0) {
                        Log (LogLevel.ERROR, "Failed for: " + i +" - "+ j);
                    }
                }
            }
            return 0;
        }


        /// <summary>
        /// Optimize the MoveSpec to maximize the CCS accuracy of the trainSet
        /// </summary>
        /// <param name="trainSet">Training examples</param>
        /// <param name="testSet">Test examples</param>
        /// <param name="algo">Viterbi or Sum-Product recursions</param>
        /// <param name="trainError">Trained error of training set</param>
        /// <param name="startParams">Start point of optimization</param>
        /// <param name="testError"></param>
        /// <returns>Trained model</returns>
        public QvModelParams OptimizeConsensus(List<CCSExample> trainSet, List<CCSExample> testSet, QvModelParams startParams, RecursionAlgo algo, int maxIterations, out float trainError, out float testError)
        {
            int iterations = 0;
            optimizer = NelderMeadOptimizer.Instance;
            Func<List<CCSExample>, Func<Vector, double>> makeObjectiveFunction = examples =>
                {
                    Func<Vector, double> f = pars =>
                        {
                            // Pipe the current spec in as the model to use in CCS.
                            var spec = ConsensusCoreWrap.QvModelParamsFromOptimizationArray(FIXED_MATCH_SCORE, pars.ToArray().Map(v => (float) v));
                            var err = LikelihoodObjective(spec, examples, algo);

                            Log(LogLevel.WARN, "Iteration {0}. Training set error: {1}", iterations, err);
                            iterations++;

                            return err;
                        };

                    return f;
                };



            // Do the full training
            var fullOptimizationStartpoint = Vector.OfEnumerable(ConsensusCoreWrap.QvModelParamsToOptimizationArray(startParams).Select(v => (double)v));
            Log("Doing full training on {0} examples", trainSet.Count);
            var fullObjectiveFunction = makeObjectiveFunction(trainSet);

            var sw = new Stopwatch();
            sw.Start();
            var result = optimizer.MinimizeFunction(fullOptimizationStartpoint, 0.005, maxIterations, fullObjectiveFunction);
            sw.Stop();

            Log("Training completed. Termination reason: {0}", result.TerminationReason);
            Log("Function evaluations: {0}", result.EvaluationCount);
            Log("Training Traces: {0}", trainSet.Count);
            Log("Iterations: {0}", iterations);
            Log("Elapsed Time: {0} mins", sw.Elapsed.TotalMinutes);
            Log("---------------------------");


            trainError = (float)result.ObjectiveValue;
            Log("Got training set error: {0}", trainError);

            var modelValues = result.X;
            Console.WriteLine(modelValues);
            var modelSpec = ConsensusCoreWrap.QvModelParamsFromOptimizationArray(FIXED_MATCH_SCORE, modelValues.ToArray().Map(v => (float) v));

            testError = 0;
            if (testSet != null)
            {
                testError = LikelihoodObjective(modelSpec, testSet, algo);
                Log("Got test set error: {0}", testError);
            }

            return modelSpec;
        }


        private CCSStream ccsAlgo = CCSStream.TrainingConfig;

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

                                var muts = GenerateMutations.GenerateMutationsForTraining(trueTpl);

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

                Log(LogLevel.WARN, "Mean LL: {0:0.0000}, Max: {1:0.0000}, Min: {2:0.0000}, Mean Aln Score: {3:0.00}", mean, maxErr, minErr, overall.Average());
                
                Log(LogLevel.INFO, "Parameters:\n{0}", 
                        Vector.OfEnumerable(ConsensusCoreWrap.QvModelParamsToArray(spec).Select(v => (double) v)).ToString());

                return mean;
            }   
        }


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
