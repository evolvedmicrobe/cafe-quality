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

            TemplatePositionTransitionGroup = new string[template.Length - 1];
            CurrentTransitionParameters = new TransitionParameters[template.Length - 1];

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
            //Fill transition probabilites appropriate for each matrix position
            FillTransitionParameters (pars);
            // clean house to be safe
            StateProbabilities.Clear();

            //We force ourselves to start and end in a match here

            // FILL THE FORWARD MATRICES
            StateProbabilities.Forward[0][0].Match = 0.0;
            StateProbabilities.Forward[0][0].Total = 0.0;
            for (int i = 1; i < (Read.Length-2); i++) {
                for (int j = 1; j < (Template.Length-2); j++) {
                    fillForwardMatrixPosition(i,j, pars);    
                }
            }

            // Now for the final transition, where we condition on ending in a match state.
            var newState = new LatentStates ();
            var transProbs = CurrentTransitionParameters [Template.Length - 1];
            var forward = StateProbabilities.Forward;
            var previous = forward [Read.Length - 2] [Template.Length - 2].Total;
            newState.Match = previous + pars.log_One_Minus_Epsilon + transProbs.log_Match;
            newState.Total = newState.Match; // Only one possible state at end
            forward [Read.Length - 1] [Template.Length - 1] = newState;


            // FILL REVERSE MATRIX
            var endi = Read.Length - 1;
            var endj = Read.Length - 1;
            StateProbabilities.Reverse[endi][endj].Match = 0.0;
            StateProbabilities.Reverse[endi][endj].Total = 0.0;
            endi--;
            endj--;
            for (int i = (endi-1); i > 0; i--) {
                for (int j = (endj-1); j > 0; j--) {
                    fillReverseMatrixPosition(i,j, pars);    
                }
            }


        }

        private void fillForwardMatrixPosition(int i, int j, ParameterSet pars)
        {
            // To store this new element
            var newState = new LatentStates ();
            var transProbs = CurrentTransitionParameters [j];

            // Match score first
            var forward = StateProbabilities.Forward;
            var previous = forward [i - 1] [j - 1].Total;
            var emissionProb = Read [i] == Template [j] ? pars.log_One_Minus_Epsilon : pars.log_Epsilon_Times_One_Third;
            newState.Match = previous + transProbs.log_Match + emissionProb;

            // Now the insertion, which is either a stick or a branch
            var probAbove = forward [i - 1] [j].Total;
            var isBranch = Template [j + 1] == Read [i];
            if (isBranch) {
                newState.Branch = probAbove + transProbs.log_Branch; // Emission probability is 1 for the same base
            } else {
                newState.Stick = probAbove + MathUtils.ONE_THIRD_LOGGED + transProbs.log_Stick;
            }

            // Now the deletion, which could be a merge or dark, in both cases the emission probability is 1
            var probLeft = forward [i] [j - 1].Total;
            // Dark first
            newState.Dark = probLeft + transProbs.log_Dark;
            if (MergePossible [j]) {
                newState.Merge = probLeft + transProbs.log_Merge; 
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

            var transProbs = CurrentTransitionParameters [j];

            // This might be a bit funky, I need to calculate for the match state the odds that we transition to every other state
            // and then combine these to get the probability.  I will be doing this by using a fake latent state and add all of them 
            var probsAfterMove = new LatentStates ();
            var reverse = StateProbabilities.Reverse;

            // The probabilities of each are the same! 

            // state -> match
            var next_match = reverse [i + 1] [j + 1].Match;
            var emissionProb = Read [i] == Template [j] ? pars.log_One_Minus_Epsilon : pars.log_Epsilon_Times_One_Third;
            probsAfterMove.Match = next_match + transProbs.log_Match + emissionProb;

            // state -> stick or state -> branch
            var isBranch = Template [j + 1] == Read [i+1];
            var next_insert = isBranch ? reverse [i + 1] [j].Branch : reverse [i + 1] [j].Stick;
            if (isBranch) {
                probsAfterMove.Branch = next_insert + transProbs.log_Branch; // Emission probability is 1 for the same base
            } else {
                probsAfterMove.Stick = next_insert + MathUtils.ONE_THIRD_LOGGED + transProbs.log_Stick;
            }

            // state -> deletion
            var next_dark = reverse [i] [j + 1].Dark;
            probsAfterMove.Dark = next_dark + transProbs.Dark;

            // state -> merge
            double next_merge;
            if (MergePossible [j]) {
                next_merge = reverse [i] [j + 1].Merge;
                probsAfterMove.Merge = transProbs.Merge + next_merge;
            }

            probsAfterMove.SetTotal ();
            probsAfterMove.Match = probsAfterMove.Total;
            probsAfterMove.Dark = probsAfterMove.Total;
//            probsAfterMove.Merge = MergePossible[j] ? probsAfterMove.Total : Double.NegativeInfinity;
//            probsAfterMove.Stick = isBranch ? Double.NegativeInfinity: probsAfterMove.Total;
//            probsAfterMove.Branch = isBranch ? probsAfterMove.Total : Double.NegativeInfinity;
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

