using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;

namespace ConstantModelOptimizer
{

    public class ReadTemplatePair
    {
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


        string Read, Template;



        public ReadTemplatePair (string read, string template)
        {
            Read = read;
            Template = template;

            if (read [0] != template [0] || read.Last () != template.Last ()) {
                throw new InvalidProgramException ("Read is expected to start in new ");
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
        public void FillMatrics(ParameterSet pars)
        {
            // Fill transition probabilites appropriate for each matrix position
            FillTransitionParameters (pars);
            // clean house to be safe
            StateProbabilities.Clear();

            //We force ourselves to start and end in a match here

            // FILL THE FORWARD MATRICES

            for (int i = 0; i < (Read.Length-1); i++) {
                for (int j = 0; j < (Template.Length-1); j++) {
                    if (j == i) {
                        Console.WriteLine ("fuck me");
                    }
                    fillForwardMatrixPosition(i,j, pars);    
                }
            }
           
            fillForwardMatrixPosition (Read.Length - 1, Template.Length - 1, pars);

            // Now for the final transition, where we condition on ending in a match state.
            var newState = new LatentStates ();
            var transProbs = CurrentTransitionParameters [Template.Length - 2];
            var forward = StateProbabilities.Forward;
            var previous = forward [Read.Length - 2] [Template.Length - 2].Total;
            newState.Match = previous + pars.log_One_Minus_Epsilon + transProbs.log_Match;
            newState.Total = newState.Match; // Only one possible state at end
            forward [Read.Length - 1] [Template.Length - 1] = newState;


            // FILL REVERSE MATRIX
            var endi = Read.Length - 2;
            var endj = Template.Length - 2;
            for (int i = endi; i >= 0; i--) {
                for (int j = endj; j >= 0; j--) {
                    fillReverseMatrixPosition(i,j, pars);    
                }
            }

            System.IO.StreamWriter sw = new System.IO.StreamWriter ("matrix.csv");
            List<Func<LatentStates, string>> grabbers = new List<Func<LatentStates, string>> () { 
                z => z.Match.ToString (),
                z => z.Stick.ToString (),
                z => z.Branch.ToString (),
                z => z.Dark.ToString (),
                z => z.Merge.ToString ()
            };
            foreach (var v in grabbers) {
                sw.WriteLine ("forward");
                var mat = StateProbabilities.Forward;
                for (int i = 0; i < mat.Length; i++) {
                    var cur = String.Join (",", mat [i].Select (x => v(x)).ToArray ());
                    sw.WriteLine (cur);
                }
                sw.WriteLine ();
//                sw.WriteLine ("reverse");
//                mat = StateProbabilities.Reverse;
//                for (int i = 0; i < mat.Length; i++) {
//                    var cur = String.Join (",", mat [i].Select (x => v(x)).ToArray ());
//                    sw.WriteLine (cur);
//                }
                sw.WriteLine ();
            }
            sw.Close();
            Console.WriteLine (StateProbabilities.Forward.Last().Last().Match + "\t" + StateProbabilities.Reverse [0] [0].Match );

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
                probsAfterMove.Dark = next_dark + transProbs.Dark;
            
                // state -> merge
                double next_merge;
                if (MergePossible [j]) {
                    next_merge = reverse [i] [j + 1].Total;
                    probsAfterMove.Merge = transProbs.Merge + next_merge;
                }
            }
            probsAfterMove.SetTotal ();
            probsAfterMove.Match = probsAfterMove.Total;
            probsAfterMove.Dark = probsAfterMove.Total;
////            probsAfterMove.Merge = MergePossible[j] ? probsAfterMove.Total : Double.NegativeInfinity;
////            probsAfterMove.Stick = isBranch ? Double.NegativeInfinity: probsAfterMove.Total;
////            probsAfterMove.Branch = isBranch ? probsAfterMove.Total : Double.NegativeInfinity;
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
            foreach (var kv in TemplatePositionTypes) {
                var noMerge = kv.Key[0]=='N';
                var cp =  pars.TransitionProbabilities[kv.Key];
                foreach (int ind in kv.Value) {
                        CurrentTransitionParameters [ind] = cp;
                }
            }           
        }
    }

}

