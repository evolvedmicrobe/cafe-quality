using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.Random;

namespace ConstantModelOptimizer
{

    public class Simulator
    {
        static Random rand = new Random();
        static char[] bases = new char[] {'A', 'G', 'C', 'T'};


        public static List<Tuple<string, string>> SimulateTemplatesAndReads(out ParameterSet pars)
        {
            List<Tuple<string, string>> pairs = new List<Tuple<string, string>> ();
            System.IO.StreamWriter sw = new System.IO.StreamWriter("Simulations.txt");
            sw.WriteLine("Template\tRead");
            pars = new ParameterSet ();
            pars.SetSingleSetDefaults ();
            for(int i=0; i < 5000; i++)
            {
                string tpl;
                string read = SimulateRead (60, pars, out tpl);
                sw.WriteLine (tpl + "\t" + read);
                pairs.Add (new Tuple<string, string> (tpl, read));
            }
            sw.Close ();
            return pairs;
        }
        /// <summary>
        /// Simulates the templates and reads.
        /// </summary>
        /// <returns>Tuple of <Template, Read> </returns>
        /// <param name="numToSimulate">Number to simulate.</param>
        /// <param name="pars">Pars.</param>
        public static List<Tuple<string, string>> SimulateTemplatesAndReads(int numToSimulate, ParameterSet pars)
        {
            List<Tuple<string, string>> pairs = new List<Tuple<string, string>> ();
            for(int i=0; i < numToSimulate; i++)
            {
                string tpl;
                string read = SimulateRead (60, pars, out tpl);
                pairs.Add (new Tuple<string, string> (tpl, read));
            }
            return pairs;
        }


        static Beta miscallRateGenerator = new Beta (1, 19);
        static Dirichlet noMergeRateGenerator = new Dirichlet(new double[] {8.0,1.0,1.0,1.0});
        static Dirichlet mergeRateGenerator = new Dirichlet(new double[] {8.0,1.0,1.0,1.0, 3.0});

        public static double SampleMisCallRate()
        {
            var rate = miscallRateGenerator.Sample ();
            return rate;
        }
        public static double[] SampleMergeRates()
        {
            var rates = mergeRateGenerator.Sample ();
            return rates;
        }
        public static double[] SampleNoMergeRates()
        {
            var rates = noMergeRateGenerator.Sample ();
            return rates;
        }
        public static string SimulateRead (int templateLength, ParameterSet pars, out string template)
        {

            template = SimulateTemplate(templateLength);

            // Now establish which parameters apply to which template positions
            var transParameters = new TransitionParameters[template.Length - 1];
            for (int j = 0; j < (template.Length - 1); j++) {
                var c1 = template [j].ToString ();
                var c2 = template [j + 1].ToString ();
                var mergePossible = c1 == c2;
                var c = mergePossible ? c1 + c2 : "N" + c2;
                if (ParameterSet.USE_DINUCLEOTIDE_MODEL) {
                    transParameters [j] = pars.TransitionProbabilities [c];
                } else {
                    transParameters [j] = mergePossible ? pars.GlobalParametersMerge : pars.GlobalParametersNoMerge;
                }
            }

            List<char> simmed = new List<char> (template.Length * 2);
            int i = 0;
            simmed.Add(template[i]);
            while (i < (template.Length - 1 )) {
                var cp = transParameters[i];
                var nextMove = SampleMultinomial(cp);
                if (nextMove == TransitionParameters.MATCH_POS)
                {
                    var nb = SampleMatchBase(template[i+1], pars.Epsilon);
                    i++;
                    simmed.Add(nb);
                }
                else if (nextMove == TransitionParameters.DARK_POS || 
                    nextMove == TransitionParameters.MERGE_POS) {
                    // have to end in a match
                    if (i < (template.Length -2)) {
                        i++;
                    }
                }
                else if(nextMove == TransitionParameters.BRANCH_POS)
                {
                    simmed.Add(template[i+1]);
                }
                else if (nextMove == TransitionParameters.STICK_POS)
                {
                    char bp;
                    var cur = template[i+1];
                    do {
                        bp = SampleBaseUniformly ();
                    } while(bp == cur);
                    simmed.Add(bp);
                }
            }
            //simmed.Add(template.Last());
            return new string(simmed.ToArray());
        }

            public static string SimulateTemplate(int templateLength = 100)
        {
            var simmed = Enumerable.Range (0, templateLength).Select (z => SampleBaseUniformly() );
            return new string (simmed.ToArray());
        }

        public static int SampleMultinomial(TransitionParameters probs)
        {
            double cdf = 0.0;
            var u = rand.NextDouble ();
            for (int i = 0; i < probs.Length; i++) {
                cdf += probs [i];
                if (u <= cdf)
                {   
                    return i;
                }
            }
            throw new Exception ("Not a valid CDF");
        }

        public static char SampleBaseUniformly()
        {
            return bases [rand.Next (4)];
        }


       

        public static char SampleMatchBase(char currentBase, double misCallProb)
        {
            var u = rand.NextDouble ();
            if (u > misCallProb) {
                return currentBase;
            }
            char toR;
            do {
                toR = SampleBaseUniformly ();
            } while(toR == currentBase);
            return toR;
        }
    }
}

