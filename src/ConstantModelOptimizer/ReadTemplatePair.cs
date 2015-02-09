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
            fillTransitionParameters (pars);
            // clean house 
            StateProbabilities.Clear();

            //We force ourselves to start and end in a math here

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
            for (int i = endi; i > 0; i--) {
                for (int j = endj; j > 0; j--) {
                    fillForwardMatrixPosition(i,j, pars);    
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
            var isBranch = Template [j + 1] = Read [i];
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
            // We want to compute the probability we were in the last state, emitted from that state, and then transitioneed here

            var newState = new LatentStates ();
            var transProbs = CurrentTransitionParameters [j];

            // Match score first
            var reverse = StateProbabilities.Reverse;
            var next = reverse [i + 1] [j + 1].Total;
            var emissionProb = Read [i] == Template [j] ? pars.log_One_Minus_Epsilon : pars.log_Epsilon_Times_One_Third;
            newState.Match = previous + transProbs.log_Match + emissionProb;


        }

        /// <summary>
        /// Fill what each transition probability is in the left-right HMM
        /// this just avoids looking up the transition parameters from a given template position each time.
        /// 
        /// A simpler solution would just use a 0 probability for each event that wasn't possible
        /// </summary>
        /// <param name="pars">Pars.</param>
        public void fillTransitionParameters(ParameterSet pars)
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

