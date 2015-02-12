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
            pars.SetDefaults ();
            foreach (var s in scorers) {
               // s.FillMatrics (pars);
            }
           
            var b = new ReadTemplatePair ("TTTT", "TTTT");
            b.FillMatrics (pars);

        }
          
    }
}
