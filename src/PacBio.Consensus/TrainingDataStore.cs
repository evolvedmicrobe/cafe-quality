using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using ConsensusCore;
using PacBio.Align;
using PacBio.IO;
using PacBio.Utils;
using PacBio.HDF;
using System.ComponentModel;
using System.Reflection;

namespace PacBio.Consensus
{
    /// <summary>
    /// A class which will hold examples of data for different partitions of the SNR/Coverage space 
    /// determined by the 
    /// </summary>
    public class TrainingDataStore
    {
        Dictionary<string, List<CCSExample>>[,] exampleStore;
        ReadConfigurationAssigner rca;
        int MaxPerReference;
        /// <summary>
        /// Indicates that we have finished adding data, and that we have shuffled all the data
        /// </summary>
        bool beenShuffledAndFullyLoaded = false;

        public TrainingDataStore (ReadConfigurationAssigner rca, int maxExamplesPerReference)
        {
            this.rca = rca;
            this.MaxPerReference = maxExamplesPerReference;
            exampleStore = new Dictionary<string, List<CCSExample>>[rca.NumberOfSNRGroups, rca.NumberOfCoverageGroups];
            for (int i = 0; i < rca.NumberOfSNRGroups; i++) {
                for (int j = 0; j < rca.NumberOfCoverageGroups; j++) {
                    exampleStore [i, j] = new Dictionary<string, List<CCSExample>> ();
                }
            }
        }

        public void PrintCounts()
        {
                Console.WriteLine(string.Join("\t", "Group","SNR","Cov","Count"));
                for (int i = 0; i < rca.NumberOfSNRGroups; i++) {
                    for (int j = 0; j < rca.NumberOfCoverageGroups; j++) {
                    try {
                        var gr = i + "-" + j;
                        var top_snr = i >= rca.MeanSNRBreakPoints.Length ? "" : " < " + rca.MeanSNRBreakPoints [i];
                        var bottom_snr = i > 0 ? rca.MeanSNRBreakPoints [i - 1] + " < " : "";
                        var snr = bottom_snr + "snr" + top_snr;

                        var top_cov = j >= rca.CoverageBreakPoints.Length ? "" : " < " + rca.CoverageBreakPoints [j];
                        var bottom_cov = j > 0 ? rca.CoverageBreakPoints [j - 1] + " < " : "";
                        var cov = bottom_cov + "cov" + top_cov;

                        var cnts = GetExamples (i, j);
                        var cnts_s = (cnts.Item1.Count + cnts.Item2.Count).ToString ();
                        Console.WriteLine(string.Join("\t", gr, snr, cov, cnts_s));
                    }
                    catch(Exception thrown) {
                        Console.WriteLine ("Print statement Failed! " + thrown.Message);
                        Console.WriteLine (thrown.StackTrace);
                    }
                }
            }
        }


        /// <summary>
        /// Add an example to the collection.
        /// </summary>
        /// <returns><c>true</c>, if example was added, <c>false</c> otherwise.</returns>
        /// <param name="example">Example.</param>
        public bool AddExample(CCSExample example)
        {
            if (beenShuffledAndFullyLoaded) {
                throw new ApplicationException ("Trying to add data after we already began taking data!");
            }
            var meanSNR = example.Trace.ZmwBases.Metrics.HQRegionSNR.Average ();
            var gr = rca.GetGroupAssignment ((float)meanSNR, example.Regions.Length);
            var result = false;
            lock (exampleStore) {
                var cur = exampleStore [gr.SnrGroup, gr.CoverageGroup];
                List<CCSExample> examples;
                bool alreadyPresent = cur.TryGetValue (example.Reference, out examples);
                if (!alreadyPresent) {
                    examples = new List<CCSExample> () { example };
                    exampleStore [gr.SnrGroup, gr.CoverageGroup] [example.Reference] = examples;
                    result = true;
                } else {
                    if (examples.Count < MaxPerReference) {
                        examples.Add (example);
                        result =  true;
                    }
                }
            }
            return result;
        }

        /// <summary>
        /// Get a training and test set for the examples specified by the SNR and Coverage grouping.
        /// </summary>
        /// <returns>The examples.</returns>
        /// <param name="snrGroup">Snr group.</param>
        /// <param name="coverageGroup">Coverage group.</param>
        public Tuple<List<CCSExample>, List<CCSExample>> GetExamples(int snrGroup, int coverageGroup, int desiredTrainingCount)
        {
            var set = GetExamples (snrGroup, coverageGroup); // Tuple<Train, Test>
            var obtainedOriginal = set.Item1.Count; 
            var obtainedBySampling = 0;
            int desiredCoverage = obtainedOriginal > 0 ? (int) Math.Round(set.Item1.Average (x => x.Regions.Length)) : -1;
            var needed = desiredTrainingCount - obtainedOriginal;
            // If we didn't get enough data, try to get "fake" data by sampling down from a higher coverage 
            // to the average of the current sample.
            // TODO: Avoid this code path
            while (set.Item1.Count < desiredTrainingCount && coverageGroup < (rca.NumberOfCoverageGroups - 2) && obtainedOriginal > 0) {
                coverageGroup++;
                var newData = GetExamples (snrGroup, coverageGroup);
                var train_n = newData.Item1.Shuffle().TakeAtMost(needed).Select(x => x.CloneWithSubSample(desiredCoverage)).ToList();
                obtainedBySampling += train_n.Count;
                var test_n = newData.Item2.Shuffle().TakeAtMost(needed).Select(x => x.CloneWithSubSample(desiredCoverage));
                set.Item1.AddRange (train_n);
                set.Item2.AddRange (test_n);
                needed = desiredTrainingCount - set.Item1.Count;
            }
            Console.WriteLine ("For " + snrGroup + " - " + coverageGroup + " From Bin: " + obtainedOriginal + " From Higher Bin: " + obtainedBySampling);
            return set;           
        }

        private Tuple<List<CCSExample>, List<CCSExample>> GetExamples(int snrGroup, int coverageGroup) {
            if (!beenShuffledAndFullyLoaded) {
                _shuffle ();
            }

            var cur = exampleStore [snrGroup, coverageGroup];
            var train = new List<CCSExample> ();
            var test = new List<CCSExample> ();
            foreach (var kv in cur) {
                var examples = kv.Value;
                var n = examples.Count / 2;
                int n_train = n + (examples.Count % 2);
                test.AddRange (examples.Take ( n_train));
                train.AddRange (examples.Skip (n_train).Take (n));
            }
            return new Tuple<List<CCSExample>, List<CCSExample>> (train, test);
        }
        private void _shuffle()
        {
            beenShuffledAndFullyLoaded = true;
            for (int i = 0; i < rca.NumberOfSNRGroups; i++) {
                for (int j = 0; j < rca.NumberOfCoverageGroups; j++) {
                    var dict = exampleStore [i, j];
                    var newDict = new Dictionary<string, List<CCSExample>> (dict.Count);
                    foreach (var kv in dict) {
                        var newL = kv.Value.Shuffle ().ToList ();
                        newDict [kv.Key] = newL;
                    }
                    exampleStore [i, j] = newDict;
                }
            }
        }
    }
}

