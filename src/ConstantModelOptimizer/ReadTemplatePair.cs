using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;

namespace ConstantModelOptimizer
{

    public class ReadTemplatePair
    {
        public string Name;

        DynamicProgrammingMatrixPair StateProbabilities;

        /// <summary>
        /// A dictionary that holds where the positions are of every dinucleotide context
        /// </summary>
        Dictionary<string, List<int>> TemplatePositionTypes;

        /// <summary>
        /// What is the dinucleotide context at position i?
        /// </summary>
        string[] TemplatePositionTransitionGroup;

        /// <summary>
        /// Arrays of ones and zeros, depending on if a merge is possible
        /// </summary>
        bool[] MergePossible;

        /// <summary>
        /// Stores the parameters for each set of transitions at the current template position, used so we don't look them up everytime
        /// </summary>
        TransitionParameters[] CurrentTransitionParameters;

        public readonly string Read; //{ get; private set; }
        public readonly string Template;// { get; private set; }


        /// <summary>
        /// This is the probability of the read at the current parameter settings.
        /// </summary>
        /// <value>The current data probability.</value>
        public double CurrentLikelihood {
            get;
            private set;
        }

        public ReadTemplatePair (string read, string template)
        {
            Read = read;
            Template = template;

            // TODO: Consider reimplementing - we want to end in a match state, 
            // which typically (for non-simulated data), means the first and last should match.
            if (read [0] != template [0] || read.Last () != template.Last ()) {
                //throw new InvalidProgramException ("Read is expected to start in new ");
            }
            // Initialize Arrays
            StateProbabilities = new DynamicProgrammingMatrixPair(read,template);
           
            // Now establish which parameters apply to which template positions
            TemplatePositionTypes = new Dictionary<string, List<int>>(ParameterSet.DiNucleotideContexts.Length);
            for (int i = 0; i < ParameterSet.DiNucleotideContexts.Length; i++) {
                var cp = ParameterSet.DiNucleotideContexts [i];
                TemplatePositionTypes [cp] = new List<int> ();
            }

            TemplatePositionTransitionGroup = new string[template.Length -1];
            CurrentTransitionParameters = new TransitionParameters[template.Length - 1];
            MergePossible = new bool[template.Length - 1 ];

            for (int i = 0; i < (template.Length - 1); i++) {
                var c1 = template [i].ToString ();
                var c2 = template [i + 1].ToString ();
                var mergePossible = c1 == c2;
                MergePossible [i] = mergePossible;

                var c = mergePossible ? c1 + c2 : "N" + c2;
                TemplatePositionTypes [c].Add (i);
                TemplatePositionTransitionGroup [i] = c;
            }
            
        }
        public double FillMatrices(ParameterSet pars)
        {
            // Fill transition probabilites appropriate for each matrix position
            FillTransitionParameters (pars);
            // clean house to be safe
            StateProbabilities.Clear();

            //We force ourselves to start and end in a match here

            // FILL THE FORWARD MATRICES
            for (int i = 0; i < (Read.Length-1); i++) {
                for (int j = 0; j < (Template.Length-1); j++) {
                    fillForwardMatrixPosition(i,j, pars);    
                }
            }
           
            fillForwardMatrixPosition (Read.Length - 1, Template.Length - 1, pars);
                  


            // FILL REVERSE MATRIX
            var endi = Read.Length - 2;
            var endj = Template.Length - 2;
            for (int i = endi; i >= 0; i--) {
                for (int j = endj; j >= 0; j--) {
                    if (i == 2 && j == 1) {
                        //    Console.WriteLine("doh!");
                    }
                    fillReverseMatrixPosition (i, j, pars);    
                }
            }
            //DumpMatrices ();


            // I don't ever save the last match value so need to add it here

            // Set the likelihood
            var likelihood = StateProbabilities.Forward.Last ().Last ().Total;
            CurrentLikelihood = likelihood;
            //Console.WriteLine (likelihood);
            // Now check for alpha beta mismatch error
            var misMatchEps = 1e-6;
            var lastMatch = Read [0] == Template [0] ? pars.log_One_Minus_Epsilon : pars.log_Epsilon_Times_One_Third;
            var beta_likelihood = (StateProbabilities.Reverse [0] [0].Total + lastMatch);
            if (Math.Abs (CurrentLikelihood - beta_likelihood) > misMatchEps) {
                throw new Exception ("Alpha-Beta mismatch error");
            }
            return CurrentLikelihood;
        }

        private void DumpMatrices()
        {
            Func<double, string> format = delegate(double x) {
                if(Double.IsNegativeInfinity(x)) {return "0";} else {return Math.Exp(x).ToString();}
                    };

            
                        System.IO.StreamWriter sw = new System.IO.StreamWriter ("matrix.csv");
                        List<Func<LatentStates, double>> grabbers = new List<Func<LatentStates, double>> () { 
                            z => z.Match,
                            z => z.Stick,
                            z => z.Branch,
                            z => z.Dark,
                            z => z.Merge, 
                            z => z.Total
                        };
                        foreach (var v in grabbers) {
                            sw.WriteLine ("forward");
                            var mat = StateProbabilities.Forward;
                            for (int i = 0; i < mat.Length; i++) {
                            var cur = String.Join (",", mat [i].Select (x => format(v(x))).ToArray ());
                                sw.WriteLine (cur);
                            }
                            sw.WriteLine ();
                            sw.WriteLine ("reverse");
                            mat = StateProbabilities.Reverse;
                            for (int i = 0; i < mat.Length; i++) {
                            var cur = String.Join (",", mat [i].Select (x => format(v(x))).ToArray ());
                                sw.WriteLine (cur);
                            }
                            sw.WriteLine ();
                        }
                        sw.Close();
        }
        private void fillForwardMatrixPosition(int i, int j, ParameterSet pars)
        {
            // To store this new element
            var newState = new LatentStates ();
            var forward = StateProbabilities.Forward;
            var matchEmissionProb = Read [i] == Template [j] ? pars.log_One_Minus_Epsilon : pars.log_Epsilon_Times_One_Third;
            var leftTransProbs = j > 0 ? CurrentTransitionParameters [j - 1] : null;
            var curTransProbs = j < CurrentTransitionParameters.Length ? CurrentTransitionParameters[j] : null;

            // Special opening case - required match
            if (i == 0 && j == 0) {
                // We are required to start in a match, so previous probability is 1
                newState.Match = matchEmissionProb;
                newState.SetTotal ();
                forward [i] [j] = newState;
                return;
            }

            // Special end case - required match 
            if ((i == Read.Length - 1) && j == (Template.Length - 1)) {
                // We are required to start in a match, so previous probability is 1
                newState.Match = forward[i-1][j-1].Total + matchEmissionProb;
                newState.SetTotal ();
                forward [i] [j] = newState;
                return;
            }

            // All others

            // Match score first
            if (i > 0 && j > 0) {
                var previous = forward [i - 1] [j - 1].Total;
                newState.Match = previous + leftTransProbs.log_Match + matchEmissionProb;
            }

            // Now the insertion, which is either a stick or a branch
            if (i > 0) {
                var probAbove = forward [i - 1] [j].Total;
                var isBranch = Template [j + 1] == Read [i];
                if (isBranch) {
                    newState.Branch = probAbove + curTransProbs.log_Branch; // Emission probability is 1 for the same base
                } else {
                    newState.Stick = probAbove + MathUtils.ONE_THIRD_LOGGED + curTransProbs.log_Stick;
                }
            }

            // Now the deletion, which could be a merge or dark, in both cases the emission probability is 1
            if (j > 0) {
                var probLeft = forward [i] [j - 1].Total;
                // Dark first
                newState.Dark = probLeft + leftTransProbs.log_Dark;
                if (MergePossible [j-1]) {
                    newState.Merge = probLeft + leftTransProbs.log_Merge; 
                }
            }
            // Now to sum it all up
            newState.SetTotal ();
            // And copy in
            forward [i] [j] = newState;

        }

        private void fillReverseMatrixPosition(int i, int j, ParameterSet pars) {
            // The backwards part of the forward backwards algorithm.  
            // We want to sum over all possible next states, the probability that we transition from 
            // this state to that state, then emit the symbol there, then all probabilities leaving that state.

            var endi = Read.Length - 2;
            var endj = Template.Length - 2;
            var matchEmissionProb = Read [i + 1] == Template [j + 1] ? pars.log_One_Minus_Epsilon : pars.log_Epsilon_Times_One_Third;
            var probsAfterMove = new LatentStates ();
            var reverse = StateProbabilities.Reverse;
            var transProbs = CurrentTransitionParameters [j];

            if (i == endi && j == endj) {
                // we have no choice but to transition and match by how the model was specified!
                probsAfterMove.Match = matchEmissionProb;
                probsAfterMove.SetTotal ();
                reverse [i] [j] = probsAfterMove;
                return;
            }

           // Now other cases

            // state -> match
            if (i < endi && j < endj) {
                var next_match = reverse [i + 1] [j + 1].Total;
                probsAfterMove.Match = next_match + transProbs.log_Match + matchEmissionProb;
            }

            // state -> stick or state -> branch
            if (i < endi) {
                var isBranch = Template [j + 1] == Read [i + 1];
                var next_insert = reverse [i + 1] [j].Total;
                if (isBranch) {
                    probsAfterMove.Branch = next_insert + transProbs.log_Branch; // Emission probability is 1 for the same base
                } else {
                    probsAfterMove.Stick = next_insert + MathUtils.ONE_THIRD_LOGGED + transProbs.log_Stick;
                }
            }

            // state -> deletion
            if (j < endj) {
                var next_dark = reverse [i] [j + 1].Total;
                probsAfterMove.Dark = next_dark + transProbs.log_Dark;
            
                // state -> merge
                double next_merge;
                if (MergePossible [j]) {
                    next_merge = reverse [i] [j + 1].Total;
                    probsAfterMove.Merge = transProbs.log_Merge + next_merge;
                }
            }
            probsAfterMove.SetTotal ();
            probsAfterMove.Match = probsAfterMove.Total;
            probsAfterMove.Dark = probsAfterMove.Total;
            probsAfterMove.Merge = probsAfterMove.Total;
            probsAfterMove.Stick = probsAfterMove.Total;
            probsAfterMove.Branch = probsAfterMove.Total;
            reverse [i] [j] = probsAfterMove;


        }

        /// <summary>
        /// Fill what each transition probability is in the left-right HMM
        /// this just avoids looking up the transition parameters from a given template position each time.
        /// 
        /// A simpler solution would just use a 0 probability for each event that wasn't possible
        /// </summary>
        /// <param name="pars">Pars.</param>
        void FillTransitionParameters(ParameterSet pars)
        {
            if (ParameterSet.USE_DINUCLEOTIDE_MODEL) {
                foreach (var kv in TemplatePositionTypes) {
                    var noMerge = kv.Key [0] == 'N';
                    var cp = pars.TransitionProbabilities [kv.Key];
                    foreach (int ind in kv.Value) {
                        CurrentTransitionParameters [ind] = cp;
                    }
                }
            } else {
                for (int i = 0; i < MergePossible.Length; i++) {
                    if (MergePossible [i]) {
                        CurrentTransitionParameters [i] = pars.GlobalParametersMerge;
                    } else {
                        CurrentTransitionParameters [i] = pars.GlobalParametersNoMerge;
                    }
                }
            }
        }

        public LatentStates GetCountsForAllColumns(ParameterSet pars, bool MergePositions) {
            if (ParameterSet.USE_DINUCLEOTIDE_MODEL) {
                throw new InvalidProgramException ("Should not be called under these conditions");
            }

            LatentStates totals = new LatentStates();
            for (int i =0; i< (Template.Length-2); i++) {
                if (MergePossible[i] == MergePositions) {
                var cnts = GetTransitionWeightsForColumn (i, pars);
                    totals.Addin (cnts);
                }
            }
            return totals;
                 
        }

        public LatentStates GetCountsForGroup(string grup, ParameterSet pars) {
            if (!ParameterSet.USE_DINUCLEOTIDE_MODEL) {
                throw new InvalidProgramException ("Should not be called under these conditions");
            }
            var columns = this.TemplatePositionTypes [grup];
            LatentStates totals = new LatentStates();
            foreach(var j in columns)
            {
                var cnts = GetTransitionWeightsForColumn (j, pars);
                totals.Addin (cnts);
            }
            return totals;
        }


        public LatentStates GetTransitionWeightsForColumn(int j, ParameterSet pars)
        {
            /* This is for the M step in the Baum-Welch algorithm. This method wants
             * to return the pseudo-counts for the expected number of times we made a K -> L
             * transition while in this column.  To do this, we will go down the entire column
             * and calculate the weighted probability we made a transition from K -> L.
             * For a move from (i,j) to (i+d, j+e), where d, l = 0 or 1 this is equal to:
             * 
             *          (F(i,j) * P_{Trans{K->L}} * Emit(i+d, j+e) * B(i+d, j+e) ) / P(Read)
             * 
             * This is sort of a mo-fo as it involves lots of progressive additions of exponentials
             */

            // Throw an exception if first or last column, these have to be matches
            if (j == Template.Length) {
                throw new Exception ("Bad column request");
            }


            var F = StateProbabilities.Forward;
            var B = StateProbabilities.Reverse;
            var transProbs = CurrentTransitionParameters [j];
            LatentStates totProbs = new LatentStates ();
            int maxRowWhereDownMovePossible = Read.Length - 3;
            int maxColumnWhereLeftMovePossible = Template.Length - 3;

            // We require the last position to be a match, so only iterate to the second to last line
            for (int i = 0; i < (Read.Length - 1); i++) {
                var curF = F [i] [j].Total;

                // Local variable to hold the weight for this particular i,j location
                var new_probs = new LatentStates ();

                // MATCH
                if (i <= maxRowWhereDownMovePossible) {
                    var emit_match = Read [i+1] == Template [j+1] ? pars.log_One_Minus_Epsilon : pars.log_Epsilon_Times_One_Third;
                    new_probs.Match = curF + transProbs.log_Match + emit_match + B [i + 1] [j + 1].Total;
                }

                // INSERTIONS, which can't account for the first or last base as that must be a match.
                // Insertions are either a stick or a branch
                if (i <= maxRowWhereDownMovePossible) {
                    var isBranch = Template [j + 1] == Read [i+1];
                    if (isBranch) {
                        new_probs.Branch = curF + transProbs.log_Branch + B[i+1][j].Total; // Emission probability is 1 for the same base
                    } else {
                        new_probs.Stick = curF + MathUtils.ONE_THIRD_LOGGED + transProbs.log_Stick + B[i+1][j].Total;
                    }
                }

                // DELETION
                // conditions already checked.
                // Dark first
                if (j <= maxColumnWhereLeftMovePossible) {
                    new_probs.Dark = curF + transProbs.log_Dark + B [i] [j + 1].Total;
                    if (MergePossible [j]) {
                        new_probs.Merge = curF + transProbs.log_Merge + B [i] [j + 1].Total; 
                    }
                }
                totProbs.Addin (new_probs);
            }
            // Now to remove P(X) from all of them (equivalent to dividing by the total probability for this path).
            totProbs.RemoveConstant (CurrentLikelihood);
            return totProbs;
        }

        public void GetInCorrectEstimate(out double numerator, out double denom) {
            // We loop from the second row to the second to last, as the first and last positions are assumed to be a match
            // In each case, F[i][j]+B[i][j] gives the total probability for a match path
            var F = StateProbabilities.Forward;
            var B = StateProbabilities.Reverse;

            numerator = Double.NegativeInfinity;
            denom = Double.NegativeInfinity;
            for (int i = 1; i < Read.Length - 1; i++) {
                for (int j = 1; j < Template.Length - 1; j++) {
                    //TODO: Move this later.
                    var amt = F [i] [j].Match + B [i] [j].Total - CurrentLikelihood;
                    if (Read [i] != Template [j]) {
                        numerator = MathUtils.logsumlog (numerator, amt);
                    }
                    denom = MathUtils.logsumlog (denom, amt);
                }
            }
//            Console.WriteLine ("Matrix Results");
//            Console.WriteLine (Math.Exp (numerator));
//            Console.WriteLine (Math.Exp (denom));
//            var sum = Enumerable.Zip (Template, Read, (z, y) => z == y ? 0 : 1).Sum();
//            Console.WriteLine (sum);
             
        }

    }

}

