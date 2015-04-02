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
            SimulateAndInfer ();
            //char bp = args [0] [0];
            //double snr = Convert.ToDouble(args[1]);
            //TrainSNRBin (bp, snr);
        }


        public static List<ReadTemplateInfo> GetSnrBinData(char bp, double snr, int numberOfSamples = 1000)
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
            var ndata = data.Shuffle ().Take(numberOfSamples).ToList();

            return ndata;

        }

        public static void TrainSNRBin(char bp, double snr) {
            Console.WriteLine ("Training " + bp.ToString () + "\t SNR: " + snr.ToString ());
            // Load 1000 samples
            var data = GetSnrBinData (bp, snr).Select(z=> new ReadTemplatePair(z.read, z.template)).ToList();
            Console.WriteLine ("Obtained " + data.Count +" samples");
            //Train on em
            var optim = new Optimizer (data);
            var pars = optim.Optimize ();
            var fname = bp.ToString () + "-" + snr.ToString () + ".txt";
            StreamWriter sw = new StreamWriter (fname);
            sw.WriteLine ("BP,SNR,Count," + pars.GetCSVHeaderLine ());
            sw.WriteLine (bp.ToString () + "," + snr.ToString () + "," +data.Count+","+ pars.GetCSVDataLine ());
            sw.Close ();
        }
        public static void InferForEdnaData()
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

        public static void SimulateAndInfer() {
            Console.WriteLine ("Hello World!");
            ParameterSet trueParameters;
            var data = Simulator.SimulateTemplatesAndReads (out trueParameters);
            System.IO.StreamWriter sw = new System.IO.StreamWriter ("TrueParametersForwardConstant.csv");
            sw.WriteLine (trueParameters.GetCSVHeaderLine ());
            sw.WriteLine (trueParameters.GetCSVDataLine ());

            //var scorers = data.Select( p=> new ReadTemplatePair(p.Item2,p.Item1)).ToList();
            // Invert the model so it's P(T|R)
            var scorers = data.Select( p=> new ReadTemplatePair(p.Item1,p.Item2)).ToList();
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

            var ll = scorers.Sum (z => z.FillMatrices (trueParameters));
            Console.WriteLine ("Real LL = " +  ll);
        }
        public class SmallRTP
        {
            public string Read, Template, Name;
        }
        public static List<SmallRTP> LoadEdnaReads()
        {
            var data = new List<SmallRTP> ();
            StreamReader sr = new StreamReader("/Users/nigel/git/cafe-quality/src/Edna/src/Edna/bin/Debug/TemplatesAndReads.txt");
            string line;
            while ((line = sr.ReadLine ()) != null) {
                string name = line.Trim();
                string template = sr.ReadLine ().Trim();
                string read = sr.ReadLine ().Trim();
                var rt = new SmallRTP () { Read = read, Template = template, Name = name };
                data.Add (rt);
            }
            return data;
        }
          
    }
}
