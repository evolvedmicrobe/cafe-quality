using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Runtime.Serialization.Formatters.Binary;
using ConsensusCore;
using PacBio.IO;
using PacBio.Utils;
using PacBio.FSharp.Utils;
using PacBio.Align;

namespace PacBio.Consensus
{
    /// <summary>
    /// Finds the consensus sequences of a group of closely related template sequences and reads drawn from them.  Assignment of reads to haplotypes is soft. That is:
    /// P(read_i is from tpl_j) = exp(Score(i,j)) / Sum_k exp(Score(i,k))
    /// where Score(i,j) is the overall Quiver alignment score of read_i on template_k.
    /// 
    /// HACKING NOTE! I'm not super happy with the clarity of this code -- there's a bit too much indexing magic that get's exposed.  
    /// Once you've grokked this code, a refactor is probably in order
    /// The basic elements are:
    ///  - The list of reads you're working with
    ///  - The set of templates you've got
    /// 
    /// Each template only get a subset of the reads mapped to it, but reads will map to multiple templates. We can calculate the posterior mapping distribution
    /// of a read against each template this way, by looking at the Quiver scores of a read on each template it is mapped to. 
    /// </summary>
    public class MultiTemplateConsensus : IDisposable
    {
        private static PacBioLogger logger = PacBioLogger.GetLogger("MultiTemplateConsensus");

        private static void Log(LogLevel level, string msg)
        {
            logger.Log(level, msg);
        }

        private static void Log(LogLevel level, string msg, params object[] args)
        {
            logger.Log(level, String.Format(msg, args));
        }

        /// <summary>
        /// Helper class to score a subset of the full set of reads against a single template
        /// </summary>
        public class SomeReadsMutationScorer : IDisposable
        {
            public SomeReadsMutationScorer(MultiReadMutationScorer initialScorer, MultiTemplateConsensus master)
            {
                Scorer = initialScorer;
                masterReadIndices = Enumerable.Range(0, initialScorer.NumReads).ToArray();
                this.Master = master;
            }

            public SomeReadsMutationScorer(SomeReadsMutationScorer parent, TrialTemplate template, int[] parentReadsToUse)
            {
                var mappedReads = parent.MappedReads;

                var newAlignedRegions = parentReadsToUse.Select(
                    idx =>
                    {
                        var originalRegAndBases = parent.Master.regions[parent.masterReadIndices[idx]];
                        var originalReg = originalRegAndBases.Item1;

                        var currentMappedRead = mappedReads[parent.masterReadIndices[idx]];

                        if (currentMappedRead == null)
                            return null;

                        var reg  = new AlignedSequenceReg(originalReg.Start, originalReg.End,
                                                          currentMappedRead.TemplateStart,
                                                          currentMappedRead.TemplateEnd, originalReg.Strand);

                        return new Tuple<AlignedSequenceReg, IZmwBases>(reg, originalRegAndBases.Item2);
                    }).Where(t => t != null).ToArray();

                Scorer = new MultiReadMutationScorer(newAlignedRegions, template, parent.Master.config);
                masterReadIndices = parentReadsToUse
                    .Select(idx => parent.masterReadIndices[idx])
                    .Where(mIdx => mappedReads[mIdx] != null)
                    .ToArray();
                this.Master = parent.Master;
            }


            /// <summary>
            /// Indexes of the reads we hold in the master read list
            /// </summary>
            private readonly int[] masterReadIndices;


            /// <summary>
            /// Master MultiTemplateConsensus caller to which this object belongs
            /// </summary>
            public MultiTemplateConsensus Master { get; private set; }

            /// <summary>
            /// The scorer for this SomeReadMutationScorer
            /// </summary>
            public MultiReadMutationScorer Scorer
            {
                get; private set;
            }
            
            public int NumReads
            {
                get { return masterReadIndices.Length; }
            }

            public TrialTemplate Template
            {
                 get { return Scorer.Template; }
            }

            public MappedRead[] MappedReads
            {
                get { return MetricSomeReadsToMaster(Scorer.MappedReads, null); }
            }

            public double[] GetBaselineScores()
            {
                return MetricSomeReadsToMaster(Scorer.GetBaselineScores(),  -Math.Sqrt(Single.MaxValue));
            }

            public int[] AllocatedEntries
            {
                get { return Scorer.AllocatedEntries; }
            }

            public int[] UsedEntries
            {
                get { return Scorer.UsedEntries; }
            }

            /// <summary>
            /// Execute Quiver template refinement on this scorer, using reads weighted according to the 
            /// </summary>
            public int ImproveConsensus(int maxIterations)
            {
                var myMappingRatios = MetricAllReadsToSome(MappingRatios());
                return FindConsensus.ImproveConsensus(Scorer, maxIterations, myMappingRatios);
            }

            /// <summary>
            /// Vector of mapping ratios for all reads, for this scorer.  Reads not mapped to this scorer at all
            /// will have a mapping ratio of 0, reads mapped to only this scorer will have a mapping ratio of Inf.
            /// </summary>
            /// <returns></returns>
            public double[] MappingRatios()
            {
                var allMappingRatios = Master.MappingRatios();
                var myIndex = Master.Scorers.FindIndex(sc => sc == this);
                return allMappingRatios[myIndex];
            }

            /// <summary>
            /// Convert an array of values computed for each read in this SomeReadsMutationScorer, 
            /// to an array over all reads in the master MultiTemplateConsensus, filling in defaultVal
            /// for reads that aren't included in this scorer.
            /// </summary>
            private T[] MetricSomeReadsToMaster<T>(T[] metric, T defaultVal)
            {
                var result = new T[Master.NumReads];

                for (int i = 0; i < result.Length; i++)
                    result[i] = defaultVal;

                for (int i = 0; i < metric.Length; i++)
                {
                    var rIdx = masterReadIndices[i];
                    result[rIdx] = metric[i];
                }

                return result;
            }

            /// <summary>
            /// Convert an array of values computed for all reads, to an array for just the reads 
            /// contained in this scorer.
            /// </summary>
            private T[] MetricAllReadsToSome<T>(T[] metric)
            {
                var result = new List<T>();

                foreach (var item in masterReadIndices)
                {
                        result.Add(metric[item]);
                }

                return result.ToArray();
            }
            
            /// <summary>
            /// Get the mutation score vector for mutation m for the reads in this scorer.
            /// The mutation scores are appropriately weighted, given the mapping ratios passed in.
            /// </summary>
            public double[] ScoreMutationWeighted(Mutation m, double[] scorerMappingRatios)
            {
                var scores = Scorer.GetScores(m); 
                var myMappingRatios = MetricAllReadsToSome(scorerMappingRatios);

                var weightedScores = scores.Map((score, i) =>
                    {
                        var r = myMappingRatios[i];

                        if (r > Math.Sqrt(float.MaxValue))
                        {
                            return scores[i];
                        }
                        else if (r < Math.Sqrt(float.Epsilon))
                        {
                            return 0.0f;
                        }
                        else
                        {
                            // Partial mapping
                            var p1 = r/(1 + r);
                            var expScore = Math.Exp(scores[i])*p1 + (1 - p1);
                            return  Math.Log(expScore);
                        }
                });

                return weightedScores;
            }

            public void Dispose()
            {
                if (Scorer != null)
                    Scorer.Dispose();
            }
        }

        // Minimum improvement in Quiver score required to split a template
        public static double splitThreshold = 6.0;

        // Minimum fraction of reads going to each side of the split
        public static double minSplitFraction = 0.1;

        // Minimum number of reads to each side of the split
        public static int minSplitReads = 8;

        // Number of haplotypes to attempt to find per recursion
        public static int maxHaplotyesPerSplit = 4;

        // Number of bases at the template ends to ignore when splitting, to help with degenerate primers
        public static int ignoreEnds = 0;

        public List<SomeReadsMutationScorer> Scorers;
        private readonly Tuple<AlignedSequenceReg, IZmwBases>[] regions;
        private readonly ScorerConfig config;

        public void Dispose()
        {
            foreach (var scorer in Scorers)
            {
                if (scorer != null)
                    scorer.Dispose();
            }
        }

        public MultiTemplateConsensus(Tuple<AlignedSequenceReg, IZmwBases>[] regions, TrialTemplate trialTemplate, ScorerConfig config)
        {
            this.config = config;
            this.regions = regions;
            
            // To start there is a single scorer which contains all the reads
            var initialScorer = new MultiReadMutationScorer(this.regions, trialTemplate, config);
            var initialSomeReads = new SomeReadsMutationScorer(initialScorer, this);
            Scorers = new List<SomeReadsMutationScorer>() { initialSomeReads };
        }

        public MultiTemplateConsensus(MultiReadMutationScorer initialScorer)
        {
            var initialSomeReads = new SomeReadsMutationScorer(initialScorer, this);
            Scorers = new List<SomeReadsMutationScorer>() { initialSomeReads };
        }

        public int NumReads
        {
            get { return regions.Length; }
        }
        

        /// <summary>
        /// Compute the mapping odds ratio for each template for each read that we hold 
        /// </summary>
        public double[][] MappingRatios()
        {
            var baselines = Scorers.Map(s => s.GetBaselineScores());

            var mappingRatios = new double[Scorers.Count][];

            for (var thisScorer = 0; thisScorer < Scorers.Count; thisScorer++)
            {
                mappingRatios[thisScorer] = new double[NumReads];

                for (var read = 0; read < NumReads; read++)
                {
                    var invMappingRatio = 0.0;

                    for (var other = 0; other < Scorers.Count; other++)
                    {
                        if (other != thisScorer)
                        {
                            invMappingRatio += Math.Exp(baselines[other][read] - baselines[thisScorer][read]);
                        }
                        
                    }
                    var ratio =  (1.0/invMappingRatio);
                    mappingRatios[thisScorer][read] = ratio;
                }
            }

            return mappingRatios;
        }

        public double[][] MappingPosteriors()
        {
            var mapRatios = MappingRatios();

            var mapPosteriors = mapRatios.Map(v => v.Map(r =>
                {                
                    if (r > Math.Sqrt(double.MaxValue))
                {
                    // Read maps essentially uniquely to this haplotype.
                    return 1.0;
                }
                    else if (r < Math.Sqrt(double.Epsilon))
                {
                    // Doesn't map strongy enough to affect outcome.
                    return 0.0;
                }
                else
                {
                    // Partial mapping
                    return r/(1.0 + r);
                }}));

            return mapPosteriors;
        }



        /// <summary>
        /// Perform quiver refinenment on the multi-template consensus
        /// </summary>
        public int[] ImproveConsensus(int maxIterations)
        {
            // Do iterGroups iterations on each template, then repeat until the total number of 
            // iteration on the longest iterating template exceeds maxIterations
            return Scorers.Select(scorer => scorer.ImproveConsensus(maxIterations)).ToArray();
        }

        /// <summary>
        /// Recursively split all the templates
        /// </summary>
        public int RecursiveSplit(int levels)
        {
            for (int l = 0; l < levels; l++)
            {
                // MultiSplit mutates Scorers, so make a local copy
                bool gotSplit = Scorers.ToList()
                    .Select(scorer => MultiSplit(scorer, maxHaplotyesPerSplit))
                    .Aggregate((a, b) => a | b);

                // Bail out if we don't split
                if (!gotSplit)
                    break;

                ImproveConsensus(1);
            }

            return Scorers.Count;
        }


        /// <summary>
        /// Split the reads in sc into multiple phases, and update Scorers accordingly.
        /// </summary>
        public bool MultiSplit(SomeReadsMutationScorer sc, int maxHaplotypes = 4)
        {
            var startPos = ignoreEnds;
            var endPos = sc.Template.Length - ignoreEnds;

            // Use all possible substitutions as potential splitting mutations
            var muts = GenerateMutations.GenerateUniqueMutations(sc.Scorer.Template)
                                        .Where(m => m.Type == MutationType.SUBSTITUTION)
                                        .Where(m => m.TemplatePosition >= startPos)
                                        .Where(m => m.TemplatePosition < endPos)
                                        .ToArray();
            
            // use weighted mutation scores - this ensures that a split primarily coming
            // from another scorer that weakly maps to this one doesn't lead to a split here.
            var mappingRatios = sc.MappingRatios();
            var mutScoreVectors = muts.Map(m => sc.ScoreMutationWeighted(m, mappingRatios));
            var positions = muts.Map(m => m.TemplatePosition);

            // Reduce the number of haplotypes to search for if there aren't enough reads.
            maxHaplotypes = (int) Math.Min(Math.Floor((float) sc.NumReads/minSplitReads), maxHaplotypes);

            var splitResult = MultiPhaseSplitter.splitHaplotypes(mutScoreVectors, positions, maxHaplotypes);

            // The splitter will return null if it can't find a good split.
            if (splitResult == null)
                return false;

            var splitScore = splitResult.Item1;
            var mutsToUse = splitResult.Item2;
            var readCounts = splitResult.Item3;
            var readFractions = splitResult.Item4;
            var readPosteriors = splitResult.Item5; 

            if (splitScore < splitThreshold)
                return false;

            var childScorers = new List<SomeReadsMutationScorer>();

            for (int i = 0; i < mutsToUse.Length; i++)
            {
                if (readCounts[i] >= minSplitReads && readFractions[i] > minSplitFraction)
                {
                    // We are going to use this split
                    // Pick all the reads the have a posterior for this 

                    // This is the subset the reads in sc that will go into the 'sub-scorer'
                    var haplotypeReadPosteriors = readPosteriors[i];

                    var readsToUse = Enumerable.Range(0, sc.NumReads).Where(idx => haplotypeReadPosteriors[idx] > 10e-8).ToArray();

                    // if the number of reads to use is less than our minimum, don't split
                    if (readsToUse.Length < minSplitReads)
                        return false;

                    var mutations = mutsToUse[i].Map(mIdx => muts[mIdx]);
                    var newTpl = sc.Template.Mutate(mutations);

                    // We are not allowed to create a new phase that matches another existing phase (other than the current sequence)
                    if (mutations.Length == 0 || sc.Master.Scorers.All(otherScorer => otherScorer.Template.Sequence != newTpl.Sequence))
                    {
                        // NOTE -- the newTpl cannot have any indels!! otherwise the template mapping of the reads will be wrong.
                        var childScorer = new SomeReadsMutationScorer(sc, newTpl, readsToUse);
                        childScorers.Add(childScorer);
                    }
                }
            }

            // If we split this set of reads into > 2 haplotypes, pull out the old haplotype and replace it with the new ones
            if (childScorers.Count > 1)
            {
                Scorers.Remove(sc);
                sc.Dispose();

                Scorers.AddRange(childScorers);

                return true;
            }

            foreach (var scorer in childScorers)
            {
                if (scorer != null)
                    scorer.Dispose();
            }

            return false;
        }


        public double[][] CoverageLevels()
        {
            var mappingPosterior = MappingPosteriors();
            var result = Scorers.Map(s => new double[s.Template.Length]);

            for (int sc = 0; sc < Scorers.Count; sc++)
            {
                var scorer = Scorers[sc];

                var reads = scorer.MappedReads;

                for (int r = 0; r < reads.Length; r++)
                {
                    var mappingWeight = mappingPosterior[sc][r];
                    var mappedRead = reads[r];

                    // This haplotype may not be using this read.
                    if (mappedRead != null)
                    {
                        var tStart = mappedRead.TemplateStart;
                        var tEnd = mappedRead.TemplateEnd;

                        for (int i = tStart; i < tEnd; i++)
                        {
                            result[sc][i] += mappingWeight;
                        }
                    }
                }
            }

            return result;
        }

        /// <summary>
        /// Compute QV and coverage statistics for each phased template sequence
        /// </summary>
		public ConsensusBaseStats[][] ConsensusStats()
		{
			var mr = MappingRatios();
			var cov = CoverageLevels();

			var result = new ConsensusBaseStats[Scorers.Count][];

			for(int i = 0; i < Scorers.Count; i++)
			{
				var sc = Scorers[i];
				var mappingRatio = mr[i];
				var coverage = cov[i];
				var tpl = sc.Template.GetSequence(Strand.Forward);

				var scores = FindConsensus.ComputeAllQVs(sc.Scorer, mappingRatio);
				var qvData = FindConsensus.ComputeConsensusQs(scores, sc.Template);

				var thisTplStats = new ConsensusBaseStats[tpl.Length];
                result[i] = thisTplStats;

				for(int p = 0; p < tpl.Length; p++)
				{
					var s = new ConsensusBaseStats
					{
						ConsensusQV = qvData.QV[p],
						Coverage = (int) Math.Round(coverage[p]),
						Base = tpl[p]
					};

					thisTplStats[p] = s;
				}
			}

            return result;
		}
    }
}
