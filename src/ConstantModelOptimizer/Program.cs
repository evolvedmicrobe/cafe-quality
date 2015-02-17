using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;


namespace ConstantModelOptimizer
{
    class MainClass
    {
        public static void Main (string[] args)
        {
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
          
    }
}
