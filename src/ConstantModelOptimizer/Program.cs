using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.IO;
using System.Diagnostics;

namespace ConstantModelOptimizer
{
    class MainClass
    {
        public static void Main (string[] args)
        {
            InferForEdnaData ();

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
            foreach (var v in data) {
                pars.SetRandomDefaults ();
                stp.Reset ();
                stp.Start ();
                var rt = new ReadTemplatePair (v.Read, v.Template);
                rt.Name = v.Name;
                var optimizer = new Optimizer(new List<ReadTemplatePair>() { rt });
                pars = optimizer.Optimize (pars);
                stp.Stop ();
                sw.WriteLine (v.Name + "," + v.Read.Length + "," + v.Template.Length + "," +stp.ElapsedMilliseconds +"," + pars.GetCSVDataLine ());
                Console.WriteLine ("Finished: " + v.Name);
                Console.WriteLine (stp.ElapsedMilliseconds.ToString ());
            }
            sw.Close ();
        }

        public static void SimulateAndInfer() {
            Console.WriteLine ("Hello World!");
            ParameterSet trueParameters;
            var data = Simulator.SimulateTemplatesAndReads (out trueParameters);
            var scorers = data.Select( p=> new ReadTemplatePair(p.Item2,p.Item1)).ToList();
            System.IO.StreamWriter sw = new System.IO.StreamWriter ("TrueParameters2.csv");
            sw.WriteLine (trueParameters.GetCSVHeaderLine ());
            sw.WriteLine (trueParameters.GetCSVDataLine ());
            sw.Close ();
            //scorers = Enumerable.Range (0, 40).Select (x => new ReadTemplatePair ("AGGT", "AGT")).ToList();
            var optim = new Optimizer (scorers);
            optim.Optimize ();

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
