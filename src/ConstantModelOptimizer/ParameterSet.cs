using System;
using System.Linq;
using System.Collections.Generic;

namespace ConstantModelOptimizer
{

    public class TransitionParametersNoHomopolymer
    {
        public double Match;
        public double Branch;
        public double Stick;
        public double Dark;

        public void Normalize()
        {
            var sum = Match+Branch+Stick+Dark;
            Match /= sum;
            Branch /= sum;
            Stick /= sum;
            Dark /= sum;
        }
    }
    public class TransitionParametersHomopolymer
    {
        public double Match;
        public double Branch;
        public double Stick;
        public double Dark;
        public double Merge;

        public void Normalize()
        {
            var sum = Match+Branch+Stick+Dark+Merge;
            Match /= sum;
            Branch /= sum;
            Stick /= sum;
            Dark /= sum;
            Merge /= sum;
        }
    }

    public class ParameterSet
    {
        /// <summary>
        /// All the di-nucleotide contexts which have their own values
        /// </summary>
        public static string[] DiNucleotideContexts = new string[] { "NA", "AA", "NG","GG", "NC", "CC", "NT", "TT"}; 


        /// <summary>
        /// DO NOT MODIFY
        /// </summary>
        public static IReadOnlyDictionary<string, int> DiNucleotideToIndex;
        /// <summary>
        /// DO NOT MODIFY
        /// </summary>
        public static IReadOnlyDictionary<int, string> IndexToDiNucleotide;

        static ParameterSet()
        {
            var diNucleotideToIndex = new Dictionary<string, int> (DiNucleotideContexts.Length);
            var indexToDiNucleotide = new Dictionary<int, string> (DiNucleotideContexts.Length);
            for(int i=0; i< ParameterSet.DiNucleotideContexts.Length; i++) {
                var c = DiNucleotideContexts [i];
                diNucleotideToIndex [c] = i;
                indexToDiNucleotide [i] = c;
            }
            DiNucleotideToIndex = diNucleotideToIndex;
            IndexToDiNucleotide = indexToDiNucleotide;
        }

        /// <summary>
        /// The probability that a match is emitted incorrectly
        /// </summary>
        public double epsilon;

        /// <summary>
        /// A dictionary which gives the probabilities of different emmissions for 
        /// </summary>
        public Dictionary<string, TransitionParametersHomopolymer> HPparameters;

        public Dictionary<string, TransitionParametersNoHomopolymer> NoHPparameters;


        public ParameterSet ()
        {


        }

        /// <summary>
        /// Set defaults, these are reasonable guesses to begin the EM step.
        /// </summary>
        public void SetDefaults()
        {
            NoHPparameters = new Dictionary<string, TransitionParametersNoHomopolymer> ();
            HPparameters = new Dictionary<string, TransitionParametersHomopolymer> ();
            // Random guesses that are roughly in line with what I think
            epsilon = 0.02;
            // I think merges happen roughly ~10% of the time based on Edna
            foreach (char c in "ACGT") {
                var s = "N" + c.ToString ();
                s = String.Intern (s);
                NoHPparameters [s] = new TransitionParametersNoHomopolymer () {
                    Match = .85,
                    Branch = 0.05,
                    Dark = 0.05,
                    Stick = 0.05
                };

                s = c.ToString () + c.ToString ();
                s = String.Intern (s);
                HPparameters [s] = new TransitionParametersHomopolymer () {
                    Match = .75,
                    Branch = 0.05,
                    Dark = 0.05,
                    Stick = 0.05,
                    Merge = 0.1
                };
            }
        }

        public string GetCSVHeaderLine()
        {
            string[] toOut = new string[29];
            toOut[0] = "Mismatch";
            var values = new string[] { "Match", "Branch", "Dark", "Stick", "Merge" };
            int j = 1;
            foreach (char c in "ACGT") {
                var s = "N" + c.ToString ();
                for (int i = 0; i < 4; i++) {
                    toOut [j] = s + "." + values [i];
                    j++;
                }
                s = c.ToString() + c.ToString ();
                for (int i = 0; i < values.Length; i++) {
                    toOut [j] = s + "." + values [i];
                    j++;
                }
            }
            return String.Join (",", toOut);
        }

        /// <summary>
        /// Gets an output line of these parameters for a CSV, matching the header output 
        /// also provided by a function call to this class.
        /// </summary>
        /// <returns>The CSV data line.</returns>
        public string GetCSVDataLine()
        {
            string[] toOut = new string[29];
            toOut[0] = "Mismatch";
            var values = new Func<TransitionParametersNoHomopolymer, double>[] { (x) => x.Match,
                x => x.Branch,
                x=>x.Dark,
                x=>x.Stick};// , "Branch", "Dark", "Stick", "Merge" };
            var values2 = new Func<TransitionParametersHomopolymer, double>[] { (x) => x.Match,
                x => x.Branch,
                x=>x.Dark,
                x=>x.Stick,
                x=>x.Merge};// , "Branch", "Dark", "Stick", "Merge" };
            int j = 1;
            foreach (char c in "ACGT") {
                var s = "N" + c.ToString ();
                var pars = NoHPparameters [s];
                for (int i = 0; i < 4; i++) {
                    toOut [j] = values [i](pars).ToString();
                    j++;
                }
                s = c.ToString() + c.ToString ();
                var pars2 = HPparameters [s];
                for (int i = 0; i < values2.Length; i++) {
                    toOut [j] = values2 [i](pars2).ToString();
                    j++;
                }
            }
            return String.Join (",", toOut);
        }
    }
}

