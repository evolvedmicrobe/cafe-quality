using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;

namespace ConstantModelOptimizer
{
    public class Simulator
    {
        static Random rand = new Random();
        static char[] bases = new char[] {'A', 'G', 'C', 'T'};


        public static void SimulateTemplatesAndReads()
        {
            System.IO.StreamWriter sw = new System.IO.StreamWriter("Simulations.txt");
            sw.WriteLine("Template\tRead");
            for(int i=0; i<20;i++)
            {
                string tpl;
                var pars = new ParameterSet ();
                pars.SetDefaults ();
                string read = SimulateRead (50, pars, out tpl);
                sw.WriteLine (tpl + "\n" + read);
                //var read = SimulateRead(50, ParameterSet.
            }
            sw.Close ();
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
                transParameters [j] = pars.TransitionProbabilities[c];
            }

            List<char> simmed = new List<char> (template.Length * 2);
            int i = 0;
            simmed.Add(template[i]);
            while (i < (template.Length - 1 )) {
                var cp = transParameters[i];
                var nextMove = SampleMultinomial(cp);
                if (nextMove == TransitionParameters.MATCH_POS)
                {
                    var nb = SampleMatchBase(template[i], pars.Epsilon);
                    i++;
                    simmed.Add(nb);
                }
                else if (nextMove == TransitionParameters.DARK_POS || 
                    nextMove == TransitionParameters.MERGE_POS) {
                    i++;
                }
                else if(nextMove == TransitionParameters.BRANCH_POS)
                {
                    simmed.Add(template[i]);
                }
                else if (nextMove == TransitionParameters.STICK_POS)
                {
                    char bp;
                    var cur = template[i];
                    do {
                        bp = SampleBaseUniformly ();
                    } while(bp != cur);
                    simmed.Add(bp);
                }
            }
            simmed.Add(template.Last());
            return new string(simmed.ToArray());
        }

        /// <summary>
        /// Duplicated code here, combine this later with the code in ReadTemplatePair.cs
        /// </summary>
        /// <param name="pars">Pars.</param>
        void FillTransitionParameters(ParameterSet pars)
        {
            foreach (var kv in TemplatePositionTypes) {
                var noMerge = kv.Key[0]=='N';
                var cp =  pars.TransitionProbabilities[kv.Key];
                foreach (int ind in kv.Value) {
                    CurrentTransitionParameters [ind] = cp;
                }
            }           
        }
        public string SimulateTemplate(int templateLength = 50)
        {
            var simmed = Enumerable.Range (0, templateLength).Select (z => SampleBaseUniformly() );
            return new string (simmed.ToArray);
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
        }

        public static char SampleBaseUniformly()
        {
            return bases [rand.Next (4)];
        }

        public static char SampleMatchBase(int currentBase, double misCallProb)
        {
            var u = rand.NextDouble ();
            if (u > misCallProb) {
                return currentBase;
            }
            char toR;
            do {
                toR = SampleBaseUniformly ();
            } while(toR != currentBase);
            return toR;
        }
    }
}

