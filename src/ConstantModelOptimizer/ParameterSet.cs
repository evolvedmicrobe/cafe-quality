﻿using System;
using System.Linq;
using System.Collections.Generic;

namespace ConstantModelOptimizer
{


    public class TransitionParameters 
    {

        public const int MATCH_POS = 0;
        public const int STICK_POS = 1;
        public const int BRANCH_POS = 2;
        public const int DARK_POS = 3;
        public const int MERGE_POS = 4;

        public double this[int i] {
            get{
                switch (i) {
                case MATCH_POS:
                    return match;
                case STICK_POS:
                    return stick;
                case BRANCH_POS:
                    return branch;
                case DARK_POS:
                    return  dark;
                case MERGE_POS:
                    return merge;
                default:
                    throw new IndexOutOfRangeException ();
                }
            }
        }

        public readonly int Length = 5;



        public double Match {
            get { return match; }
            set {
                match = value;
                log_match = Math.Log (value);
            }
        }
        public double Branch {
            get { return branch; }
            set {
                branch = value;
                log_branch = Math.Log (value);
            }
        }
        public double Stick
        {
            get { return stick; }
            set {
                stick = value;
                log_stick = Math.Log (value);
            }
        }
        public double Dark
        {
            get { return dark; }
            set {
                dark = value;
                log_dark = Math.Log (value);
            }
        }
        /// <summary>
        /// Zero if not in a homopolymer context
        /// </summary>
        public double Merge {
            get { return merge; }
            set {
                merge = value;
                log_merge = Math.Log (value);
            }
        }


        public double log_Match {
            get { return log_match; }
        }
        public double log_Branch {
            get { return log_branch; }

        }
        public double log_Stick
        {
            get { return log_stick; }

        }
        public double log_Dark
        {
            get { return log_dark; }

        }
        /// <summary>
        /// Zero if not in a homopolymer context
        /// </summary>
        public double log_Merge {
            get { return log_merge; }

        }


        private double log_match, log_branch, log_stick, log_dark, log_merge;
        private double match, branch, stick, dark, merge;

        public TransitionParameters MakeCopy()
        {
            var ntp = new TransitionParameters ();
            ntp.Merge = Merge;
            ntp.Stick = Stick;
            ntp.Match = Match;
            ntp.Branch = Branch;
            ntp.Dark = Dark;
            return ntp;
        }
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
        #region STATICS
        /// <summary>
        /// Whether to have one set of parameters for each dinucleotide context or not.
        /// </summary>
        public const bool USE_DINUCLEOTIDE_MODEL = true;


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
        #endregion
        /// <summary>
        /// The probability that a match is emitted incorrectly
        /// </summary>
        public double Epsilon
        {
            get { return epsilon; }
            set {
                log_one_minus_epsilon = Math.Log (1-value);
                log_epsilon_times_one_third = Math.Log (value) + MathUtils.ONE_THIRD_LOGGED;
                epsilon = value;
            }
        }

        public double log_Epsilon_Times_One_Third {
            get { return log_epsilon_times_one_third;}
        }

        public double log_One_Minus_Epsilon {
            get{ return log_one_minus_epsilon; }
        }

        private double epsilon;
        private double log_one_minus_epsilon;
        private double log_epsilon_times_one_third;

        /// <summary>
        /// The parameters to use if we only have one set, 
        /// and are not dividing by di-nucleotide context.
        /// </summary>
        public TransitionParameters GlobalParametersNoMerge;
        public TransitionParameters GlobalParametersMerge;



        /// <summary>
        /// A dictionary which gives the probabilities of different emmissions for 
        /// </summary>
        public Dictionary<string, TransitionParameters> TransitionProbabilities;

        /// <summary>
        /// Makes a new parameter set by deep copying another.
        /// </summary>
        /// <param name="toCopy">To copy.</param>
        public ParameterSet(ParameterSet toCopy)
        {
            this.Epsilon = toCopy.Epsilon;
            TransitionProbabilities = new Dictionary<string, TransitionParameters> ();
            foreach (var kv in toCopy.TransitionProbabilities) {
                this.TransitionProbabilities [kv.Key] = kv.Value.MakeCopy ();
            }
        }

        public ParameterSet ()
        {


        }

        /// <summary>
        /// Set defaults, these are reasonable guesses to begin the EM step.
        /// </summary>
        public void SetDefaults()
        {
            if (USE_DINUCLEOTIDE_MODEL) {
                TransitionProbabilities = new Dictionary<string, TransitionParameters> ();
                // Random guesses that are roughly in line with what I think
                Epsilon = 0.04;
                // I think merges happen roughly ~10% of the time based on Edna
                foreach (char c in "ACGT") {
                    var s = "N" + c.ToString ();
                    s = String.Intern (s);
                    TransitionProbabilities [s] = new TransitionParameters () {
                        Match = .95,
                        Branch = 0.05,
                        Dark = 0.05,
                        Stick = 0.05
                    };
                    TransitionProbabilities [s].Normalize ();
                    s = c.ToString () + c.ToString ();
                    s = String.Intern (s);
                    TransitionProbabilities [s] = new TransitionParameters () {
                        Match = .75,
                        Branch = 0.05,
                        Dark = 0.05,
                        Stick = 0.05,
                        Merge = 0.1
                    };
                    TransitionProbabilities [s].Normalize ();
                }
            } else {
                Epsilon = 0.04;
                GlobalParametersNoMerge = new TransitionParameters () {
                    Match = 1.0,
                    Branch = 0,
                    Dark = 0,
                    Stick = 0
                };
                GlobalParametersNoMerge.Normalize ();
                GlobalParametersMerge = new TransitionParameters () {
                    Match = 1.0,
                    Branch = 0.0,
                    Dark = 0.0,
                    Stick = 0.0,
                    Merge = 0.0
                };
                GlobalParametersMerge.Normalize ();
            }
        }

        public void SetRandomDefaults()
        {
            Epsilon = Simulator.SampleMisCallRate ();
            if (USE_DINUCLEOTIDE_MODEL) {
                TransitionProbabilities = new Dictionary<string, TransitionParameters> ();
                // Random guesses that are roughly in line with what I think
               
                // I think merges happen roughly ~10% of the time based on Edna
                foreach (char c in "ACGT") {
                    var s = "N" + c.ToString ();
                    s = String.Intern (s);
                    var np = Simulator.SampleNoMergeRates ();

                    TransitionProbabilities [s] = new TransitionParameters () {
                        Match = np[0],
                        Branch = np[1],
                        Dark = np[2],
                        Stick = np[3]
                    };
                    TransitionProbabilities [s].Normalize ();
                    s = c.ToString () + c.ToString ();
                    s = String.Intern (s);
                    np = Simulator.SampleMergeRates ();
                    TransitionProbabilities [s] = new TransitionParameters () {
                        Match = np[0],
                        Branch = np[1],
                        Dark = np[2],
                        Stick = np[3],
                        Merge = np[4]
                    };
                    TransitionProbabilities [s].Normalize ();
                }
            } else {
                var np = Simulator.SampleNoMergeRates ();
                GlobalParametersNoMerge = new TransitionParameters () {
                    Match = np[0],
                    Branch = np[1],
                    Dark = np[2],
                    Stick = np[3]
                };
                GlobalParametersNoMerge.Normalize ();
                np = Simulator.SampleMergeRates ();
                GlobalParametersMerge = new TransitionParameters () {
                    Match = np[0],
                    Branch = np[1],
                    Dark = np[2],
                    Stick = np[3],
                    Merge = np[4]
                };
                GlobalParametersMerge.Normalize ();
            }
        }

        public string GetCSVHeaderLine()
        {
            string[] toOut = new string[41];
            toOut[0] = "Mismatch";
            var values = new string[] { "Match", "Branch", "Dark", "Stick", "Merge" };
            int j = 1;
            foreach (char c in "ACGT") {
                var s = "N" + c.ToString ();
                for (int i = 0; i < values.Length; i++) {
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
            string[] toOut = new string[41];
            toOut[0] = Epsilon.ToString();
            var values2 = new Func<TransitionParameters, double>[] { 
                x => x.Match,
                x => x.Branch,
                x=>x.Dark,
                x=>x.Stick,
                x=>x.Merge};// , "Branch", "Dark", "Stick", "Merge" };
            int j = 1;
            foreach (char c in "ACGT") {
                var s = "N" + c.ToString ();
                var pars = TransitionProbabilities [s];
                for (int i = 0; i < values2.Length; i++) {
                    toOut [j] = values2 [i](pars).ToString();
                    j++;
                }
                s = c.ToString() + c.ToString ();
                var pars2 = TransitionProbabilities [s];
                for (int i = 0; i < values2.Length; i++) {
                    toOut [j] = values2 [i](pars2).ToString();
                    j++;
                }
            }
            return String.Join (",", toOut);
        }
    }
}
