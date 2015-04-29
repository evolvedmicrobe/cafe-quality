using System;
using NUnit.Framework;
using ConsensusCore;
using ConstantModelOptimizer;
//using PacBio.Consensus;
using System.Linq;

namespace Tests
{
    using TransitionParameters = ConstantModelOptimizer.TransitionParameters;

    [TestFixture ()]
    public class SnrTests
    {
       
        [Test()]
        public void TestWeirdRead() {
            var read = "CTGGGGAT";
            var tpl =  "CTGGGGGAT";
            var StateProbabilities = new DynamicProgrammingMatrixPair(read,tpl);
            TransitionParameters[] tps = new TransitionParameters[]{ 
                new TransitionParameters(0.9661447, 0.006205411, 0.019176495, 0.008473387),
                new TransitionParameters(0.8382683, 0.045386146, 0.009378846, 0.106966725),
                new TransitionParameters(0.7079446, 0.010019695, 0.014392839, 0.267642900),
                new TransitionParameters(0.7023541, 0.008000416, 0.014005718, 0.275639804),
                new TransitionParameters(0.8072699, 0.021285437, 0.016855127, 0.154589508),
                new TransitionParameters(0.8733340, 0.010229981, 0.014839555, 0.101596418),
                new TransitionParameters(0.9580190, 0.012116631, 0.022381979, 0.007482362),
                new TransitionParameters(0.9657584, 0.008621670, 0.021174486, 0.004445403),
                new TransitionParameters(0.9445740, 0.009206869, 0.038266323, 0.007952813)
            };
            ReadTemplatePair rtp = new ReadTemplatePair ( read, tpl);
            ParameterSet ps = new ParameterSet ();
            ps.Epsilon = 0.002552835;
            rtp.FillMatrices (tps, ps);
            var incorrect = rtp.CurrentLikelihood;
            Console.WriteLine (rtp.CurrentLikelihood);
            rtp = new ReadTemplatePair (tpl, tpl);
            rtp.FillMatrices (tps, ps);
            var correct = rtp.CurrentLikelihood;
            Console.WriteLine (incorrect - correct);


        }

        [Test()]
        public void TestMultiReadMutationScorer()
        {
            //First let's get some parameters
            SNR snr = new SNR (10.0, 7.0, 5.0, 11.0);
            ParameterSet ps = new ParameterSet ();
            ps.Epsilon = 0.002671256; //HARD CODED IN C++
            foreach (var ctx in ParameterSet.DiNucleotideContexts) {
                var pars = ContextParameterProvider.GetTransitionParameters(ctx, snr);
                var n_params = new ConstantModelOptimizer.TransitionParameters ();
                n_params.Match = Math.Exp (pars.Match);
                n_params.Dark = Math.Exp (pars.Deletion);
                n_params.Branch = Math.Exp (pars.Branch);
                n_params.Stick = Math.Exp (pars.Stick);
                n_params.Merge = 0;
                ps.TransitionProbabilities [ctx] = n_params;
            }

            var read = "ACGTACGT";
            var template = "ACGTCGT";

            ReadTemplatePair rtp = new ReadTemplatePair (read, template);
            rtp.FillMatrices (ps);
            //rtp.DumpMatrices ("test1.csv");
            //4.9422203069979727
            var res = rtp.CurrentLikelihood;
            Console.WriteLine (res);

            ReadTemplatePair rtpTemp = new ReadTemplatePair ("ACGTACGT", "ACGTACGT");
            rtpTemp.FillMatrices (ps);
            rtpTemp.DumpMatrices ("test1.csv");
            //-0.584415070238446
            var resTemp = rtpTemp.CurrentLikelihood;
            Console.WriteLine (resTemp);

            rtpTemp = new ReadTemplatePair ("ACGTACGT", "ACGTGCGT");
            rtpTemp.FillMatrices (ps);
            rtpTemp.DumpMatrices ("test1.csv");
            resTemp = rtpTemp.CurrentLikelihood;
            Console.WriteLine (resTemp);


           

            template = read;
            rtp = new ReadTemplatePair (read, template);
            rtp.FillMatrices (ps);
            // C# wants this to be -0.584415070238446
            res = rtp.CurrentLikelihood;
            Console.WriteLine (res);
            //
            template = "ACCTCGT";
            rtp = new ReadTemplatePair (read, template);
            rtp.FillMatrices (ps);
            //rtp.DumpMatrices ();
            // C# wants this to be -8.49879694901693
            var res2 = rtp.CurrentLikelihood;


            template = "ACGTGT";
            rtp = new ReadTemplatePair (read, template);
            rtp.FillMatrices (ps);
            var res4 = rtp.CurrentLikelihood;



            template = "AACGTCGT";
            rtp = new ReadTemplatePair (read, template);
            rtp.FillMatrices (ps);
            var res5 = rtp.CurrentLikelihood;



            //rtp.DumpMatrices ();
            Console.WriteLine (res2);

        }

        [Test()]
        public void CompareToCPlusPlus()
        {
            //First let's get some parameters
            SNR snr = new SNR (10.0, 7.0, 5.0, 11.0);
            ParameterSet ps = new ParameterSet ();
            ps.Epsilon = 0.002671256; //HARD CODED IN C++
            foreach (var ctx in ParameterSet.DiNucleotideContexts) {
                var pars = ContextParameterProvider.GetTransitionParameters(ctx, snr);
                var n_params = new ConstantModelOptimizer.TransitionParameters ();
                n_params.Match = Math.Exp (pars.Match);
                n_params.Dark = Math.Exp (pars.Deletion);
                n_params.Branch = Math.Exp (pars.Branch);
                n_params.Stick = Math.Exp (pars.Stick);
                n_params.Merge = 0;
                ps.TransitionProbabilities [ctx] = n_params;
            }

            var read = "ACGTACGT";
            var template = "ACGTCGT";

            ReadTemplatePair rtp = new ReadTemplatePair (read, template);
            rtp.FillMatrices (ps);
            //rtp.DumpMatrices ("testF.csv");
            var res = rtp.CurrentLikelihood;
            Console.WriteLine (res);

            template = read;
            rtp = new ReadTemplatePair (read, template);
            rtp.FillMatrices (ps);
            // C# wants this to be -0.584415070238446
            res = rtp.CurrentLikelihood;
            Console.WriteLine (res);


            rtp = new ReadTemplatePair ("ACCTCGT","ACGTCGT");
            rtp.FillMatrices (ps);
            //rtp.DumpMatrices ();
            res = rtp.CurrentLikelihood;
            Console.WriteLine (res);
//
            template = "ACCTCGT";
            rtp = new ReadTemplatePair (read, template);
            rtp.FillMatrices (ps);
            //rtp.DumpMatrices ();
            // C# wants this to be -8.49879694901693
            var res2 = rtp.CurrentLikelihood;


            template = "ACGTGT";
            rtp = new ReadTemplatePair (read, template);
            rtp.FillMatrices (ps);
            var res4 = rtp.CurrentLikelihood;



            template = "AACGTCGT";
            rtp = new ReadTemplatePair (read, template);
            rtp.FillMatrices (ps);
            var res5 = rtp.CurrentLikelihood;



            //rtp.DumpMatrices ();
            Console.WriteLine (res2);

        }

        [Test ()]
        public void TestSnrCalibration ()
        {                
            //SNR settings = new SNR (6.0, 6.0, 6.0, 6.0);
            SNR settings = new SNR (10.0, 7.0, 5.0, 11.0);
            var result2 = ContextParameterProvider.GetTransitionParameters ("NC", settings);

            // Compare to values calculated in R.
            var result = ContextParameterProvider.GetTransitionParameters ("NA", settings);
            double eps = 0.001;
            Assert.That (Math.Abs(result.Branch - -3.2866373) < eps);
            Assert.That (Math.Abs (-2.7102493 - result.Deletion) < eps);
            Assert.That (Math.Abs (-0.1375533 - result.Match) < eps);
            Assert.That (Math.Abs (-3.7044989 - result.Stick) < eps);
            settings.Dispose ();
            //Now let's check to make sure the right channel is used.
            settings = new SNR (0.0, 8.0, 0.0, 0.0);
            result = ContextParameterProvider.GetTransitionParameters ("NC", settings);
            Assert.That (Math.Abs(result.Branch - -3.63717488) < eps);
            Assert.That (Math.Abs (-3.55678567 - result.Deletion) < eps);
            Assert.That (Math.Abs (-0.09390532 - result.Match) < eps);
            Assert.That (Math.Abs (-3.35888384 - result.Stick) < eps);
            settings.Dispose ();

        }


//        public void WriteAlphaBeta(SparseQvSumProductMultiReadMutationScorer scorer) {
//            var data = scorer.Read (0);
//            System.IO.StreamWriter writer = new System.IO.StreamWriter ("Matrices.csv");
//            foreach (var alpha in new bool[] {true, false}) {
//
//                var mat = alpha ? scorer.AlphaMatrix (0) : scorer.BetaMatrix (0);
//                var nrow = mat.Rows ();
//                var ncol = mat.Columns ();
//
//                for (int x = 0; x < nrow; x++) {
//                    var row = Enumerable.Range (0, ncol)
//                    .Select (y => mat.IsAllocated (x, y) ? mat.Get (x, y) : float.NaN);
//                    writer.WriteLine (row.Select (v => v.ToString ()));
//                }
//            }
//            writer.Close ();
//        }
//

        [Test()]
        public void TestMutationScorer() {
           // MatrixTester mt = new MatrixTester ();
            //var res = mt.TestMutationScorer ();
            //Assert.IsTrue (res == 0);
        }

        [Test()]
        public void TestLikelihoodCalculation()
        {
            //First let's get some parameters
            SNR snr = new SNR (10.0, 7.0, 5.0, 11.0);
            ParameterSet ps = new ParameterSet ();
            ps.Epsilon = 0.002671256; //HARD CODED IN C++
            foreach (var ctx in ParameterSet.DiNucleotideContexts) {
                var pars = ContextParameterProvider.GetTransitionParameters(ctx, snr);
                var n_params = new ConstantModelOptimizer.TransitionParameters ();
                n_params.Match = Math.Exp (pars.Match);
                n_params.Dark = Math.Exp (pars.Deletion);
                n_params.Branch = Math.Exp (pars.Branch);
                n_params.Stick = Math.Exp (pars.Stick);
                ps.TransitionProbabilities [ctx] = n_params;
            }

            // Now to simulate some data
            var data = Simulator.SimulateTemplatesAndReads (1, ps).First();
            var tpl = data.Item1;
            var read_str = data.Item2;

            // Now get the consensus core score 
            var ctx_params = new ContextParameters (snr);
            var scoreDiff = 18;
            var fastScoreThreshold = -12.5;
            var bo = new BandingOptions ( scoreDiff);
            var qc = new QuiverConfig(ctx_params, bo, fastScoreThreshold);

            double cplusplusScore = 1;

            SparseQvSumProductMultiReadMutationScorer scorer = new SparseQvSumProductMultiReadMutationScorer (qc, tpl);
            using (var read = new Read ("Test", read_str)) {
                using (var mappedRead = new MappedRead (read, StrandEnum.FORWARD_STRAND, 0, tpl.Length,
                                        true, true)) { //TODO: Figure out about this add read suff

                    var active = scorer.AddRead (mappedRead);
                    //WriteAlphaBeta (scorer);
                    Assert.IsTrue (active == AddReadResult.SUCCESS);
                    cplusplusScore = scorer.BaselineScore();
                }
            }
            // And the C# score to compare.
            ReadTemplatePair rtp = new ReadTemplatePair (read_str, tpl);
            rtp.FillMatrices (ps);
            double csharpScore = rtp.CurrentLikelihood;
            Assert.AreEqual (csharpScore, cplusplusScore);

        }
    }
}

