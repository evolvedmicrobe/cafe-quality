using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;

namespace ConstantModelOptimizer
{
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
        string[] TemplatePositionMapping;

        public ReadTemplatePair (string read, string template)
        {
            if (read [0] != template [0] || read.Last () != template.Last ()) {
                throw new InvalidProgramException ("Read is expected to start in new ");
            }
            // Initialize Arrays
            Merge = new DynamicProgrammingMatrixPair(read,template);
            Branch = new DynamicProgrammingMatrixPair(read, template);
            Stick = new DynamicProgrammingMatrixPair(read, template);
            Dark = new DynamicProgrammingMatrixPair(read, template);
            Match = new DynamicProgrammingMatrixPair(read, template);

            // Now establish which parameters apply to which template positions
            TemplatePositionTypes = new Dictionary<string, List<int>>(ParameterSet.DiNucleotideContexts.Length);
            for (int i = 0; i < ParameterSet.DiNucleotideContexts.Length; i++) {
                var cp = ParameterSet.DiNucleotideContexts [i];
                TemplatePositionTypes [cp] = new List<int> ();
            }
            TemplatePositionMapping = new string[template.Length - 1];
            for (int i = 0; i < (template.Length - 1); i++) {
                var c1 = template [i].ToString ();
                var c2 = template [i + 1].ToString ();
                var c = c1 == c2 ? c1 + c2 : "N" + c2;
                TemplatePositionTypes [c].Add (i);
                TemplatePositionMapping [i] = c;
            }
        }
        public void FillMatrics(ParameterSet pars)
        {


        }
    }

}

