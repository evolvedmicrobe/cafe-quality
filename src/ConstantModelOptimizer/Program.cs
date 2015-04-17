using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.IO;
using System.Diagnostics;
using PacBio.Utils;

namespace ConstantModelOptimizer
{
    class MainClass
    {
        public static void Main (string[] args)
        {
            
            // Simulate and infer
            SimulateAndInfer ();
            
            // TRAIN A SNR BIN
            //char bp = args [0] [0];
            //double snr = Convert.ToDouble(args[1]);
            //TrainSNRBin (bp, snr);
            // return;

            // Test a weird read
            TestWeirdRead ();
            return;

            //var model = "ACGTAA";
            //var outcome = "ACGTAA";
            //var pars = new ParameterSet ();
            //pars.SetUniform (0.05);
            //var pair = new ReadTemplatePair (outcome, model);
            //var ll = pair.FillMatrices (pars);
            //double n, d;
            //pair.GetInCorrectEstimate (out n, out d);
            //Console.WriteLine ("ll = {0}, n = {1}, d = {2}", ll, Math.Exp (n), Math.Exp (d));
            //for (int i = 0; i < model.Length - 1; i += 1) {
            //    Console.WriteLine (pair.GetTransitionWeightsForColumn (i, pars));
            //}
            //return;
            
            //SimulateAndInfer (Int32.Parse (args [0]), Int32.Parse (args [1]));
            //return;

        }


        public static void TestWeirdRead ()
        {
            var mdl = "CTGGGGGAT";
            //var StateProbabilities = new DynamicProgrammingMatrixPair(mdl, o1);
            TransitionParameters[] tps = new TransitionParameters[] { 
                new TransitionParameters (0.9661447, 0.006205411, 0.019176495, 0.008473387),
                new TransitionParameters (0.8382683, 0.045386146, 0.009378846, 0.106966725),
                new TransitionParameters (0.7079446, 0.010019695, 0.014392839, 0.267642900),
                new TransitionParameters (0.7023541, 0.008000416, 0.014005718, 0.275639804),
                new TransitionParameters (0.8072699, 0.021285437, 0.016855127, 0.154589508),
                new TransitionParameters (0.8733340, 0.010229981, 0.014839555, 0.101596418),
                new TransitionParameters (0.9580190, 0.012116631, 0.022381979, 0.007482362),
                new TransitionParameters (0.9657584, 0.008621670, 0.021174486, 0.004445403),
                // new TransitionParameters(0.9445740, 0.009206869, 0.038266323, 0.007952813)
            };
            string o1 = "CTGGGGAT", o2 = "CTGGGGGAT";
            ReadTemplatePair rtp = new ReadTemplatePair (o1, mdl);
            ParameterSet ps = new ParameterSet ();
            ps.Epsilon = 0.002552835;
            rtp.FillMatrices (tps, ps);
            var ll1 = rtp.CurrentLikelihood;
            Console.WriteLine ("{0}: {1}", o1, ll1);
            rtp.DumpMatrices (Console.Out);
            rtp = new ReadTemplatePair (o2, mdl);
            rtp.FillMatrices (tps, ps);
            var ll2 = rtp.CurrentLikelihood;
            Console.WriteLine ("{0}: {1}", o2, ll2);
            rtp.DumpMatrices (Console.Out);
        }


        public static List<ReadTemplateInfo> GetSnrBinData (char bp, double snr, int numberOfSamples = 1000)
        {

            // Define function to test if SNR is in range.
            Func<double, bool> inRange = (x) => Math.Abs (x - snr) < 0.5;
            Func<ZmwInfo, bool> selectSNR;
            switch (bp) {
            case 'A':
                selectSNR = (x) => inRange (x.SnrA);
                break;
            case 'G':
                selectSNR = x => inRange (x.SnrG);
                break;
            case 'C':
                selectSNR = x => inRange (x.SnrC);
                break;
            case 'T':
                selectSNR = x => inRange (x.SnrT);
                break;
            default:
                throw new ApplicationException ();
            }

            // Get a list of all holes in that range.
            var snrData = LoadData.LoadSNRs ();
            HashSet<int> acceptableHoles = new HashSet<int> (snrData.Where (z => selectSNR (z)).Select (p => p.HoleNumber));

            // Now load all such reads
            var data = LoadData.LoadSampleData ().Where (z => acceptableHoles.Contains (z.Hole)).ToList ();
            var ndata = data.Shuffle ().Take (numberOfSamples).ToList ();

            return ndata;

        }

        public static void TrainSNRBin (char bp, double snr)
        {
            Console.WriteLine ("Training " + bp.ToString () + "\t SNR: " + snr.ToString ());
            // Load 1000 samples
            var data = GetSnrBinData (bp, snr).Select (z => new ReadTemplatePair (z.read, z.template)).ToList ();
            Console.WriteLine ("Obtained " + data.Count + " samples");
            //Train on em
            var optim = new Optimizer (data);
            var pars = optim.Optimize ();
            var fname = bp.ToString () + "-" + snr.ToString () + ".txt";
            StreamWriter sw = new StreamWriter (fname);
            sw.WriteLine ("BP,SNR,Count," + pars.GetCSVHeaderLine ());
            sw.WriteLine (bp.ToString () + "," + snr.ToString () + "," + data.Count + "," + pars.GetCSVDataLine ());
            sw.Close ();
        }

        public static void InferForEdnaData ()
        {
            Console.WriteLine ("Hello World!");
            var data = LoadEdnaReads ();
            var pars = new ParameterSet ();
            pars.SetRandomDefaults ();
            StreamWriter sw = new StreamWriter ("EdnaFit.csv");
            sw.WriteLine ("Name, ReadLength, TemplateLength, ComputeMS," + pars.GetCSVHeaderLine ());
            Stopwatch stp = new Stopwatch ();
            int out_count = 0;
            foreach (var v in data) {
                List<ReadTemplatePair> splits;
                if (v.Template.Length > 200) {
                    splits = AlignmentSplitter.SplitSequence (v.Read, v.Template);
                } else {
                    splits = new List<ReadTemplatePair> () { new ReadTemplatePair (v.Read, v.Template) };
                }
                if (splits != null) {
                    pars.SetRandomDefaults ();
                    stp.Reset ();
                    stp.Start ();
                    //var rt = new ReadTemplatePair (v.Read, v.Template);
                    //rt.Name = v.Name;
                    var optimizer = new Optimizer (splits);
                    pars = optimizer.Optimize (pars);
                    stp.Stop ();
                    sw.WriteLine (v.Name + "," + v.Read.Length + "," + v.Template.Length + "," + stp.ElapsedMilliseconds + "," + pars.GetCSVDataLine ());
                    sw.Flush ();
                    out_count++;
                    Console.WriteLine ("Finished: " + v.Name);
                    Console.WriteLine (stp.ElapsedMilliseconds.ToString ());
                } else {
                    Console.WriteLine ("Malformed read/template pair");
                }
            }
            Console.WriteLine (out_count);
            sw.Close ();
        }

        public static void SimulateAndInfer (int l, int n)
        {
            Console.WriteLine ("Hello World!");
            ParameterSet trueParameters;
            var data = Simulator.SimulateTemplatesAndReads (out trueParameters, l);
            var scorers = data.Select (p => new ReadTemplatePair (p.Item2, p.Item1)).ToList ();
            var ll = scorers.Sum (z => z.FillMatrices (trueParameters));
            System.IO.StreamWriter sw = new System.IO.StreamWriter ("TrueParameters2.csv");
            sw.WriteLine ("Likelihood," + trueParameters.GetCSVHeaderLine ());
            sw.WriteLine (ll.ToString () + "," + trueParameters.GetCSVDataLine ());
            sw.Close ();
            //scorers = Enumerable.Range (0, 40).Select (x => new ReadTemplatePair ("AGGT", "AGT")).ToList();

            var optim = new Optimizer (scorers);
            var fit = optim.Optimize ();
            sw.WriteLine (fit.GetCSVDataLine());

            // Now try it the other way, P(R|T)
            scorers = data.Select( p=> new ReadTemplatePair(p.Item2,p.Item1)).ToList();
            optim = new Optimizer (scorers);
            fit = optim.Optimize ();
            sw.WriteLine (fit.GetCSVDataLine());
            sw.Close ();

            ll = scorers.Sum (z => z.FillMatrices (trueParameters));
            Console.WriteLine ("Real LL = " + ll);
            System.IO.StreamWriter dw = new System.IO.StreamWriter ("SimData.cpp");
            dw.WriteLine ();
            dw.WriteLine ("#include <iostream>");
            dw.WriteLine ("#include <vector>");
            dw.WriteLine ();
            dw.WriteLine ("#include <boost/lexical_cast.hpp>");
            dw.WriteLine ();
            dw.WriteLine ("#include <unitem/Fit.hpp>");
            dw.WriteLine ("#include <unitem/Types.hpp>");
            dw.WriteLine ();
            dw.WriteLine ();
            dw.WriteLine ("int main(int argc, char **argv)");
            dw.WriteLine ("{");
            dw.WriteLine ("    using namespace std;");
            dw.WriteLine ();
            dw.WriteLine ("    const size_t trunc = (argc > 1) ? boost::lexical_cast<size_t>(argv[1]) : 1;");
            dw.WriteLine ("    const double substitution = {0};", trueParameters.Epsilon);
            dw.WriteLine ();
            int m = 1;
            foreach (var p in scorers.TakeAtMost(n)) {
                double pll = p.FillMatrices (trueParameters);
                p.DumpMatrices (String.Format ("matrices{0}.csv", m));
                dw.WriteLine ("    {");
                dw.WriteLine ("        const Outcome o{0}(\"{1}\");", m, p.Read);
                dw.WriteLine ("        vector<PrTransD> m{0};", m);
                for (int i = 0; i < p.CurrentTransitionParameters.Length; i++) {
                    var tp = p.CurrentTransitionParameters [i];
                    dw.WriteLine ("        m{0}.push_back(PrTransD({1}, {2}, {3}, {4}, {5}));", m, p.Template [i], tp.Match, tp.Branch, tp.Stick, tp.Dark);
                }
                dw.WriteLine ("        m{0}.push_back(PrTransD({1}, 0.25, 0.25, 0.25, 0.25));", m, p.Template [p.CurrentTransitionParameters.Length]);
                dw.WriteLine ("        const double csll{0} = {1};", m, pll);
                dw.WriteLine ("        cerr << \"fitting model {0} ..\" << endl;", m);
                dw.WriteLine ("        auto r{0} = FitOutcomeToModel(o{0}, ModelD(m{0}), substitution, trunc);", m);
                //dw.WriteLine ("    cout << \"csLL = \" << csll{0} << \", cppLL = \" << r{0}.LgLikelihood << endl;", m);
                dw.WriteLine ("        if (trunc == 1 && fabs(csll{0} - r{0}.LgLikelihood) > ConstantsD::Eps)", m);
                dw.WriteLine ("            cerr << \"ll mismatch! csLL = \" << csll{0} << \", cppLL = \" << r{0}.LgLikelihood << endl;", m);
                dw.WriteLine ("    }");
                dw.WriteLine ();
                m++;
            }
            dw.WriteLine ("    return 0;");
            dw.WriteLine ("}");
            dw.Close ();
        }

        public class SmallRTP
        {
            public string Read, Template, Name;
        }

        public static List<SmallRTP> LoadEdnaReads ()
        {
            var data = new List<SmallRTP> ();
            StreamReader sr = new StreamReader ("/Users/nigel/git/cafe-quality/src/Edna/src/Edna/bin/Debug/TemplatesAndReads.txt");
            string line;
            while ((line = sr.ReadLine ()) != null) {
                string name = line.Trim ();
                string template = sr.ReadLine ().Trim ();
                string read = sr.ReadLine ().Trim ();
                var rt = new SmallRTP () { Read = read, Template = template, Name = name };
                data.Add (rt);
            }
            return data;
        }
          
    }
}
