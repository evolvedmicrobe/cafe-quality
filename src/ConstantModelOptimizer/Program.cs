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
            var data = Simulator.SimulateTemplatesAndReads ();
            var scorers = data.Select( p=> new ReadTemplatePair(p.Item2,p.Item1)).ToList();
            var pars = new ParameterSet ();
            //scorers = Enumerable.Range (0, 40).Select (x => new ReadTemplatePair ("AGGT", "AGT")).ToList();
            pars.SetDefaults ();
            var optim = new Optimizer (scorers);
            optim.Optimize ();
            foreach (var s in scorers) {
//                s.FillMatrices (pars);
            }
           
            var b = new ReadTemplatePair ("ACGT", "AGT");
            b.FillMatrices (pars);

        }
          
    }
}
