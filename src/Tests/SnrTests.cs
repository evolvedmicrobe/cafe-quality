using System;
using NUnit.Framework;
using ConsensusCore;
using ConstantModelOptimizer;
using PacBio.Consensus;
using System.Linq;

namespace Tests
{
    [TestFixture ()]
    public class SnrTests
    {
       

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
            //rtp.DumpMatrices ();
            var res = rtp.CurrentLikelihood;
            Console.WriteLine (res);

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
            rtp.DumpMatrices ();
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



            rtp.DumpMatrices ();
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
            //rtp.DumpMatrices ();
            var res = rtp.CurrentLikelihood;
            Console.WriteLine (res);

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
            rtp.DumpMatrices ();
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



            rtp.DumpMatrices ();
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


        public void WriteAlphaBeta(SparseQvSumProductMultiReadMutationScorer scorer) {
            var data = scorer.Read (0);
            System.IO.StreamWriter writer = new System.IO.StreamWriter ("Matrices.csv");
            foreach (var alpha in new bool[] {true, false}) {

                var mat = alpha ? scorer.AlphaMatrix (0) : scorer.BetaMatrix (0);
                var nrow = mat.Rows ();
                var ncol = mat.Columns ();

                for (int x = 0; x < nrow; x++) {
                    var row = Enumerable.Range (0, ncol)
                    .Select (y => mat.IsAllocated (x, y) ? mat.Get (x, y) : float.NaN);
                    writer.WriteLine (row.Select (v => v.ToString ()));
                }
            }
            writer.Close ();
        }


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
            var diag_cross = 4;
            var scoreDiff = 18;
            var fastScoreThreshold = -12.5;
            var bo = new BandingOptions (diag_cross, scoreDiff);
            var qc = new QuiverConfig(ctx_params, bo, fastScoreThreshold);

            double cplusplusScore = 1;

            SparseQvSumProductMultiReadMutationScorer scorer = new SparseQvSumProductMultiReadMutationScorer (qc, tpl);
            using (var read = new Read ("Test", read_str)) {
                using (var mappedRead = new MappedRead (read, StrandEnum.FORWARD_STRAND, 0, tpl.Length,
                                        true, true)) { //TODO: Figure out about this add read suff

                    var active = scorer.AddRead (mappedRead);
                    WriteAlphaBeta (scorer);
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

