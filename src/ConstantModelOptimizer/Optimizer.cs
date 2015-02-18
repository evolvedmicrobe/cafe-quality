using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;


namespace ConstantModelOptimizer
{
    public class Optimizer
    {

        public List<ReadTemplatePair> data;
        public Optimizer (List<ReadTemplatePair> data)
        {
            this.data = data;
        }
        public ParameterSet Optimize(ParameterSet startParms = null)
        {
            ParameterSet pars;
            if (startParms == null) {
                pars = new ParameterSet ();
                pars.SetRandomDefaults ();
            } else {
                pars = new ParameterSet (startParms);
            }


            double ll = double.MinValue;
            double ll_dif = double.MaxValue;
            double termination_dif = 1e-2;
            System.IO.StreamWriter sw = new System.IO.StreamWriter ("Parameters2.csv");
            sw.WriteLine ("Likelihood," + pars.GetCSVHeaderLine ());
            bool first = false;
            while (ll_dif > termination_dif) {

                // Fill forward-backward matrices and get likelihood
                // E - Step

                Parallel.ForEach( data, z => z.FillMatrices(pars));
                var new_ll = data.Sum(z => z.CurrentLikelihood);
                // Output the initial parameters
                if (first) {
                    sw.WriteLine (new_ll + "," + pars.GetCSVDataLine ());
                    first = false;
                }
                //var new_ll = data.Sum(z => z.FillMatrices(pars));
                Console.WriteLine ("Log likelihood: \t" + new_ll);
                if (new_ll < ll) {
                    // In EM the likelihood always goes up!
                    // I have observed some cases where the likelihoood only changes by ~1e-4%, 
                    // and I think this is due to numerical issues with pseudo counts near the optimum.
                    // Hence, the second condition listed above.
                    if ((1.0 - new_ll / ll) > 1e-4) {
                        break;
                    }
                    throw new ApplicationException ("Someone didn't code the algorithm correctly");
                }
                ll_dif = new_ll - ll;
                ll = new_ll;
                // Update the parameters
                // M - Step

                // Update transition probabilities
                Console.WriteLine ("ctx\tMatch\tStick\tBranch\tDark\tMerge");
                if (ParameterSet.USE_DINUCLEOTIDE_MODEL) {
                    Parallel.ForEach (ParameterSet.DiNucleotideContexts, ctx => { 
//                    foreach (var ctx in ParameterSet.DiNucleotideContexts) {
                        // Get all the pseudo-counts across datasets

                        var cnts = new LatentStates ();
                        foreach (var d in data) {
                            cnts.Addin (d.GetCountsForGroup (ctx, pars));
                        }
                        // Now divide by the total to get the expected transition types
                        cnts.SetTotal ();
                        cnts.RemoveConstant (cnts.Total);
                        //if (ctx == "AG") {

                        var cp = pars.TransitionProbabilities [ctx];
                        cp.Match = Math.Exp (cnts.Match);
                        cp.Stick = Math.Exp (cnts.Stick);
                        cp.Branch = Math.Exp (cnts.Branch);
                        cp.Dark = Math.Exp (cnts.Dark);
                        cp.Merge = Math.Exp (cnts.Merge);
                        Console.WriteLine (String.Join ("\t", ctx, cp.Match.ToString (), cp.Stick.ToString (), cp.Branch.ToString (), cp.Dark.ToString (), cp.Merge.ToString ()));
                        //}
                    }
                    );
                } else {
                    var cnts = new LatentStates ();
                    foreach (var d in data) {
                        cnts.Addin(d.GetCountsForAllColumns(pars,true));
                    }
                    cnts.SetTotal ();
                    cnts.RemoveConstant (cnts.Total);
                    var cp = pars.GlobalParametersMerge;
                    cp.Match = Math.Exp (cnts.Match);
                    cp.Stick = Math.Exp (cnts.Stick);
                    cp.Branch = Math.Exp (cnts.Branch);
                    cp.Dark = Math.Exp (cnts.Dark);
                    cp.Merge = Math.Exp (cnts.Merge);
                    Console.WriteLine (String.Join ("\t", "Merge", cp.Match.ToString (), cp.Stick.ToString (), cp.Branch.ToString (), cp.Dark.ToString (), cp.Merge.ToString ()));

                    cnts = new LatentStates ();
                    foreach (var d in data) {
                        cnts.Addin(d.GetCountsForAllColumns(pars,false));
                    }
                    cnts.SetTotal ();
                    cnts.RemoveConstant (cnts.Total);
                    cp = pars.GlobalParametersNoMerge;
                    cp.Match = Math.Exp (cnts.Match);
                    cp.Stick = Math.Exp (cnts.Stick);
                    cp.Branch = Math.Exp (cnts.Branch);
                    cp.Dark = Math.Exp (cnts.Dark);
                    cp.Merge = Math.Exp (cnts.Merge);
                    Console.WriteLine (String.Join ("\t", "No-Merge", cp.Match.ToString (), cp.Stick.ToString (), cp.Branch.ToString (), cp.Dark.ToString (), cp.Merge.ToString ()));
                }
                // Update the miscall probability
                double pseudo_incorrect = Double.NegativeInfinity; 
                double pseudo_total = Double.NegativeInfinity;
                foreach (var d in data) {
                    double num, denom;
                    d.GetInCorrectEstimate (out num, out denom);
                    pseudo_incorrect = MathUtils.logsumlog (pseudo_incorrect, num);
                    pseudo_total = MathUtils.logsumlog (pseudo_total, denom);
                    //Console.WriteLine (Math.Exp (num) / Math.Exp(denom));
                    //correct = MathUtils.logsumlog (correct, d.GetCorrectEstimate ());
                }


                var incorrect = Math.Exp(pseudo_incorrect - pseudo_total);
                Console.WriteLine("Eps\t" + incorrect);
                pars.Epsilon = incorrect;
                if (ParameterSet.USE_DINUCLEOTIDE_MODEL) {
                    sw.WriteLine (new_ll + "," +  pars.GetCSVDataLine ());
                }
                sw.Flush ();
            }
            sw.Close ();
            return pars;
        }

    }
}

