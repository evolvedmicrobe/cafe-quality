using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;

namespace ConstantModelOptimizer
{
    public struct TransitionParameters {
        public TransitionParametersNoHomopolymer NoHPparams;
        public TransitionParametersHomopolymer HPparams;
    }

    public class ReadTemplatePair
    {
        DynamicProgrammingMatrixPair Stick;
        DynamicProgrammingMatrixPair Branch;
        DynamicProgrammingMatrixPair Merge;
        DynamicProgrammingMatrixPair Dark;
        DynamicProgrammingMatrixPair Match;

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

        List<DynamicProgrammingMatrixPair> Matrices = new List<DynamicProgrammingMatrixPair>();

        string Read, Template;

        public void 


        public ReadTemplatePair (string read, string template)
        {
            Read = read;
            Template = template;

            if (read [0] != template [0] || read.Last () != template.Last ()) {
                throw new InvalidProgramException ("Read is expected to start in new ");
            }
            // Initialize Arrays
            Merge = new DynamicProgrammingMatrixPair(read,template);
            Branch = new DynamicProgrammingMatrixPair(read, template);
            Stick = new DynamicProgrammingMatrixPair(read, template);
            Dark = new DynamicProgrammingMatrixPair(read, template);
            Match = new DynamicProgrammingMatrixPair(read, template);

            Matrices = new List<DynamicProgrammingMatrixPair> () { Merge, Branch, Stick, Dark, Match };

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
            //Fill matrix
            fillTransitionParameters (pars);
            // clean house 
            foreach (var mat in Matrices) {
                mat.Clear ();
            }

            //We force ourselves to start and end in a math here
            Match.forward [0] [0] = 0.0;

            // Now let's fill the forward matrices
            for (int i = 1; i < Read.Length; i++) {
                for (int j = 1; j < Template.Length; j++) {
                            
                }
            }

        }

        private void fillMatrix(int i, int j)
        {
            // Match score first!
            if 

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
                var cp = noMerge ? pars.NoHPparameters [kv.Key] : null;
                var cp2 = noMerge ? null : pars.HPparameters [kv.Key];
                foreach (int ind in kv.Value) {
                    if (noMerge) {
                        CurrentTransitionParameters [ind].NoHPparams = cp;
                    } else {
                        CurrentTransitionParameters [ind].HPparams = cp2;
                    }
                }
            }
           
        }
    }

}

