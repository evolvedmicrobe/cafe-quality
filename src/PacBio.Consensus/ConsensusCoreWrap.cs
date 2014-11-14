using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.IO.Compression;
using System.Linq;
using System.Runtime.Serialization.Formatters.Binary;
using ConsensusCore;
using PacBio.Align;
using PacBio.IO;
using PacBio.Utils;

namespace PacBio.Consensus
{

    public class ConsensusCorePoa
    {
        public static string FindConsensus(string[] sequences, out float consensusScore)
        {
            var config = new PoaConfig(false);

            // Put reads into POA
            var stringVect = new StringVector(sequences);

            // Compute consensus
            var p = PoaConsensus.FindConsensus(stringVect, config);

            consensusScore = p.Score();
            return p.Sequence();
        }
    }
    
    /// <summary>
    /// Utility code for getting data into ConsensusCore
    /// </summary>
    public class ConsensusCoreWrap
    {
        public static QvSequenceFeatures MakeQvSequenceFeatures(int startBase, int length, IZmwBases b)
        {
            using (FloatArray delQv = MakeFloatArray(b.DeletionQV, startBase, length),
                              delTag = MakeFloatArray(b.DeletionTag, startBase, length),
                              insQv = MakeFloatArray(b.InsertionQV, startBase, length),
                              subsQv = MakeFloatArray(b.SubstitutionQV, startBase, length),
                              mergeQv = MakeFloatArray(b.MergeQV, startBase, length)
                )
            {
                return new QvSequenceFeatures(b.Sequence.Substring(startBase, length), insQv.cast(), subsQv.cast(),
                                              delQv.cast(), delTag.cast(), mergeQv.cast());
            }
        }

        public static QvSequenceFeatures MakeQvSequenceFeatures(int referenceChunkStart, int referenceChunkSize,
                                                                IAlnSummary alnSummary)
        {
            var alignment = alnSummary.ReadAlignment();

            var useCell =
                alignment.Cells.Select(
                    c =>
                    c.Template >= referenceChunkStart && c.Template < referenceChunkStart + referenceChunkSize &&
                    c.Mode != AlignMode.Delete).ToArray();


            using (FloatArray delQv = MakeFloatArray(FilterArray(alnSummary.ReadMetric<byte>("DeletionQV"), useCell)),
                              delTag = MakeFloatArray(FilterArray(alnSummary.ReadMetric<sbyte>("DeletionTag"), useCell)),
                              insQv = MakeFloatArray(FilterArray(alnSummary.ReadMetric<byte>("InsertionQV"), useCell)),
                              subsQv = MakeFloatArray(FilterArray(alnSummary.ReadMetric<byte>("SubstitutionQV"), useCell)),
                              mergeQv = MakeFloatArray(FilterArray(alnSummary.ReadMetric<byte>("MergeQV"), useCell)))
            {
                var sequence = new string(FilterArray(alignment.ReadMarkup().ToCharArray(), useCell));
                Debug.Assert(!(sequence.Contains("-")));

                return new QvSequenceFeatures(sequence, insQv.cast(), subsQv.cast(),
                                              delQv.cast(), delTag.cast(), mergeQv.cast());
            }
        }

        private static T[] FilterArray<T>(IEnumerable<T> array, IList<bool> filter)
        {
            return array.Where((v, idx) => filter[idx]).ToArray();
        }



        private static FloatArray MakeFloatArray(IList<float> a, int start = 0, int length = -1)
        {
            if (length < 0)
                length = a.Count;

            var fa = new FloatArray(length);

            for (var i = 0; i < length; i++)
                fa.setitem(i, a[start + i]);

            return fa;

        }


        private static FloatArray MakeFloatArray(IList<byte> a, int start = 0, int length = -1)
        {
            if (length < 0)
                length = a.Count;

            var fa = new FloatArray(length);
            for (int i = 0; i < length; i++)
                fa.setitem(i, a[start + i]);

            return fa;
        }


        private static FloatArray MakeFloatArray(IList<char> a, int start = 0, int length = -1)
        {
            if (length < 0)
                length = a.Count;

            var fa = new FloatArray(length);

            for (int i = 0; i < length; i++)
                fa.setitem(i, a[start + i]);

            return fa;
        }


        private static FloatArray MakeFloatArray(IList<sbyte> a, int start = 0, int length = -1)
        {
            if (length < 0)
                length = a.Count;

            var fa = new FloatArray(length);
            for (var i = 0; i < length; i++)
                fa.setitem(i, a[start + i]);

            return fa;
        }


        public static float[] QvModelParamsToArray(QvModelParams p)
        {
            var mergeArray = FloatArray.frompointer(p.Merge);
            var mergeSArray = FloatArray.frompointer(p.MergeS);

            var d = new[]
            {
                p.Match,
                p.Mismatch,
                p.MismatchS,
                p.Branch,
                p.BranchS,
                p.DeletionN,
                p.DeletionWithTag,
                p.DeletionWithTagS,
                p.Nce,
                p.NceS,
                mergeArray.getitem(0),
                mergeArray.getitem(1),
                mergeArray.getitem(2),
                mergeArray.getitem(3),
                mergeSArray.getitem(0),
                mergeSArray.getitem(1),
                mergeSArray.getitem(2),
                mergeSArray.getitem(3)
            };

            return d;
        }


        public static QvModelParams QvModelParamsFromArray(float[] arr)
        {
            var m = new QvModelParams(
                arr[0],
                arr[1],
                arr[2],
                arr[3],
                arr[4],
                arr[5],
                arr[6],
                arr[7],
                arr[8],
                arr[9],
                arr[10],
                arr[11],
                arr[12],
                arr[13],
                arr[14],
                arr[15],
                arr[16],
                arr[17]
            );

            return m;
        }
    }

    public class ScorerConfig : IDisposable
    {
        public QuiverConfigTable Parameters;
        public RecursionAlgo Algorithm = RecursionAlgo.Viterbi;
        public bool HasChemistryOverride = false;

        public void Dispose()
        {
            if (Parameters != null)
                Parameters.Dispose();
        }
    }


    /// <summary>
    /// Computes the scores of a mutations with respect to an initial template for a set of reads.  This is essentially a wrapper around the MultiReadMutationScorer
    /// class in ConsensusCore.  Some things, notably the Mutation type and the 'TrialTemplate' type have one type on the C# side and a different type in ConsensusCore/SWIG.
    /// There may be some opportunity to clean that up.
    /// </summary>
    public class MultiReadMutationScorer : IDisposable
    {
        private static PacBioLogger logger = PacBioLogger.GetLogger("MultiReadMutationScorer");

        private static void Log(LogLevel level, string msg)
        {
            logger.Log(level, msg);
        }

        private static void Log(LogLevel level, string msg, params object[] args)
        {
            logger.Log(level, String.Format(msg, args));
        }

        public float ScoreDiff = 15;
        public float DynamicAdjustFactor = 0.0f;
        public float DynamicAdjustOffset = -0.3f;
        #if DIAGNOSTIC
        public const float AddThreshold = 0.15f;
        #else
        public const float AddThreshold = 1.0f;
        #endif

        /// <summary>
        /// Number of bases of the adapter included at front of sequence
        /// </summary>
        private int StartAdapterBases;

        /// <summary>
        /// Number of bases of adapter included at the end of sequence
        /// </summary>
        private int EndAdapterBases;

        private AbstractMultiReadMutationScorer scorer;

        /// <summary>
        /// Construct a mutation evalutor for the pulse observations in encapsulated by read, on template TrialTemplate.
        /// </summary>
        public MultiReadMutationScorer(IEnumerable<AlignedSequenceReg> regions, IZmwBases bases,
                                       TrialTemplate trialTemplate, ScorerConfig config)
        {
            StartAdapterBases = trialTemplate.StartAdapterBases;
            EndAdapterBases = trialTemplate.EndAdapterBases;

            var strandTpl = trialTemplate.GetSequence(Strand.Forward);

            if (config.Algorithm == RecursionAlgo.Viterbi)
            {
                scorer = new SparseSseQvMultiReadMutationScorer(config.Parameters, strandTpl);
            }
            else if (config.Algorithm == RecursionAlgo.Prob)
            {
                scorer = new SparseSseQvSumProductMultiReadMutationScorer(config.Parameters, strandTpl);
            }
            else
            {
                throw new ApplicationException("unrecognized recursion algorthim");
            }

            foreach (var r in regions)
            {
                var name = TraceReference.CreateSpringfieldSubread(bases, r.Start, r.End);
                var chem = config.HasChemistryOverride ? "*" : bases.Zmw.Movie.SequencingChemistry;

                using (var qsf = ConsensusCoreWrap.MakeQvSequenceFeatures(r.Start, r.End - r.Start, bases))
                using (var read = new Read(qsf, name, chem))
                using (var mappedRead = new MappedRead(read, (StrandEnum) r.Strand, r.TemplateStart, r.TemplateEnd,
                                                       r.AdapterHitBefore, r.AdapterHitAfter))
                {
                    if (!scorer.AddRead(mappedRead, AddThreshold))
                    {
                        #if DIAGNOSTIC
                        scorer.AddRead(mappedRead, 1.0f);

                        var prefix = name.Replace("/", ".");

                        WriteAlphaBetaMatrices(readBitmap.Count, prefix);
                        #else
                        Log(LogLevel.DEBUG, "Skipped adding read '{0}' -- more than {1}% of entries used.", name, AddThreshold * 100);
                        #endif
                    }
                }
            }
        }


        /// <summary>
        /// Construct a mutation evalutor for the pulse observations in encapsulated by read, on template TrialTemplate.
        /// </summary>
        public MultiReadMutationScorer(IEnumerable<Tuple<AlignedSequenceReg, IZmwBases>> regions,
                                       TrialTemplate trialTemplate, ScorerConfig config)
        {
            StartAdapterBases = trialTemplate.StartAdapterBases;
            EndAdapterBases = trialTemplate.EndAdapterBases;

            var strandTpl = trialTemplate.GetSequence(Strand.Forward);

            if (config.Algorithm == RecursionAlgo.Viterbi)
            {
                scorer = new SparseSseQvMultiReadMutationScorer(config.Parameters, strandTpl);
            }
            else if (config.Algorithm == RecursionAlgo.Prob)
            {
                scorer = new SparseSseQvSumProductMultiReadMutationScorer(config.Parameters, strandTpl);
            }
            else
            {
                throw new ApplicationException("unrecognized recursion algorthim");
            }

            foreach (var regAndBases in regions)
            {
                var r = regAndBases.Item1;
                var bases = regAndBases.Item2;
                var name = TraceReference.CreateSpringfieldSubread(bases, r.Start, r.End);
                var chem = config.HasChemistryOverride ? "*" : bases.Zmw.Movie.SequencingChemistry;

                using (var qsf = ConsensusCoreWrap.MakeQvSequenceFeatures(r.Start, r.End - r.Start, bases))
                using (var read = new Read(qsf, name, chem))
                using (var mappedRead = new MappedRead(read, (StrandEnum)r.Strand, r.TemplateStart, r.TemplateEnd,
                                            r.AdapterHitBefore, r.AdapterHitAfter))
                {
                    if (!scorer.AddRead(mappedRead, AddThreshold))
                    {
                        #if DIAGNOSTIC
                        scorer.AddRead(mappedRead, 1.0f);

                        var prefix = name.Replace("/", ".");

                        WriteAlphaBetaMatrices(readBitmap.Count, prefix);
                        #else
                        Log(LogLevel.DEBUG, "Skipped adding read '{0}' -- more than {1}% of entries used.", name, AddThreshold * 100);
                        #endif
                    }
                }
            }

        }


        /// <summary>
        /// Construct a mutation evalutor for the pulse observations in encapsulated by read, on template TrialTemplate.
        /// </summary>
        public MultiReadMutationScorer(IEnumerable<IAlnSummary> alignments, int referenceChunkStart,
                                       TrialTemplate trialTemplate, ScorerConfig config)
        {
            StartAdapterBases = trialTemplate.StartAdapterBases;
            EndAdapterBases = trialTemplate.EndAdapterBases;

            var strandTpl = trialTemplate.GetSequence(Strand.Forward);

            if (config.Algorithm == RecursionAlgo.Viterbi)
            {
                scorer = new SparseSseQvMultiReadMutationScorer(config.Parameters, strandTpl);
            }
            else if (config.Algorithm == RecursionAlgo.Prob)
            {
                scorer = new SparseSseQvSumProductMultiReadMutationScorer(config.Parameters, strandTpl);
            }
            else
            {
                throw new ApplicationException("unrecognized recursion algorthim");
            }

            var referenceChunkSize = trialTemplate.Sequence.Length;
            var referenceChunkEnd = referenceChunkStart + referenceChunkSize;

            foreach (var a in alignments)
            {
                // we have to subtract out the current referenceChunk start because the provided
                // trialTemplate has already been chunked.
                // TODO (lhepler) -- fix how this is handled, it is very inelegant at this moment
                var templateStart = Math.Max(a.TemplateStart, referenceChunkStart) - referenceChunkStart;
                var templateEnd = Math.Min(a.TemplateEnd, referenceChunkEnd) - referenceChunkStart;
                var chem = config.HasChemistryOverride ? "*" : a.SequencingChemistry;

                using (var qsf = ConsensusCoreWrap.MakeQvSequenceFeatures(referenceChunkStart, referenceChunkSize, a))
                using (var read = new Read(qsf, a.ReadName, chem))
                using (var mappedRead = new MappedRead(read, (StrandEnum)a.Strand, templateStart, templateEnd))
                {
                    if (!scorer.AddRead(mappedRead, AddThreshold))
                    {
                        #if DIAGNOSTIC
                        scorer.AddRead(mappedRead, 1.0f);

                        var prefix = name.Replace("/", ".");

                        WriteAlphaBetaMatrices(readBitmap.Count, prefix);
                        #else
                        Log(LogLevel.DEBUG, "Skipped adding read '{0}' -- more than {1}% of entries used.", a.ReadName, AddThreshold * 100);
                        #endif
                    }
                }
            }
        }


        /// <summary>
        /// Gets the number of allocated cells for each added MappedRead
        /// </summary>
        public int[] AllocatedEntries
        {
            get
            {
                using (var allocatedMatrixEntries = scorer.AllocatedMatrixEntries())
                {
                    var r = new int[allocatedMatrixEntries.Count];
                    allocatedMatrixEntries.CopyTo(r);
                    return r;
                }
            }
        }

        /// <summary>
        /// Gets the used cells for each added MappedRead
        /// </summary>
        public int[] UsedEntries
        {
            get
            {
                using (var usedMatrixEntries = scorer.UsedMatrixEntries())
                {
                    var r = new int[usedMatrixEntries.Count];
                    usedMatrixEntries.CopyTo(r);
                    return r;
                }
            }
        }


        private void WriteCsvMatrix(int i, bool alpha, string prefix)
        {
            using (var stream = new FileStream(prefix + ".csv.gz", FileMode.Create))
            using (var gzstream = new GZipStream(stream, CompressionMode.Compress))
            using (var writer = new StreamWriter(gzstream))
            {
                var mat = alpha ? scorer.AlphaMatrix(i) : scorer.BetaMatrix(i);
                var nrow = mat.Rows();
                var ncol = mat.Columns();

                for (int x = 0; x < nrow; x++)
                {
                    var row = Enumerable.Range(0, ncol)
                        .Select(y => mat.IsAllocated(x, y) ? mat.Get(x, y) : float.NaN);
                    writer.WriteLine(row.Select(v => v.ToString()));
                }
            }
        }


        public void WriteAlphaBetaMatrices(int i, string prefix)
        {
            WriteCsvMatrix(i, true, prefix + ".alpha");
            WriteCsvMatrix(i, false, prefix + ".beta");
        }


        /// <summary>
        /// Estimate of the memory consumption in MB of the recursion matrices
        /// </summary>
        public float AllocatedSizeMB
        {
            get
            {
                return  AllocatedEntries.Sum() * 4.0f / 1e6f;
            }
        }

        /// <summary>
        /// Number of forward-backward flip-flops taken to converge for each read
        /// </summary>
        public int[] NumFlipFlops
        {
            get
            {
                using (var nFlipFlops = scorer.NumFlipFlops())
                {
                    var r = new int[nFlipFlops.Count];
                    nFlipFlops.CopyTo(r);
                    return r;
                }
            }
        }

        /// <summary>
        /// Access the ConsensusCore MappedRead object for each read
        /// </summary>
        public MappedRead[] MappedReads
        {
            get
            {
                var r = new MappedRead[scorer.NumReads()];

                for (int i = 0; i < scorer.NumReads(); i++)
                {
                    r[i] = scorer.Read(i);
                }

                return r;
            }
        }

        public int NumReads
        {
            get { return scorer.NumReads(); }
        }

        /// <summary>
        /// Copy constructor
        /// </summary>
        private MultiReadMutationScorer(MultiReadMutationScorer other)
        {
            StartAdapterBases = other.StartAdapterBases;
            EndAdapterBases = other.EndAdapterBases;

            if (other.scorer is SparseSseQvMultiReadMutationScorer)
            {
                scorer = new SparseSseQvMultiReadMutationScorer((SparseSseQvMultiReadMutationScorer) other.scorer);
            }
            else if (other.scorer is SparseSseQvSumProductMultiReadMutationScorer)
            {
                scorer = new SparseSseQvSumProductMultiReadMutationScorer((SparseSseQvSumProductMultiReadMutationScorer) other.scorer);
            }
            else
            {
                throw new ApplicationException("Unrecognized MultiReadMutationScorer type");
            }
        }

        /// <summary>
        /// Make a clone of the scorer with the same set of reads.
        /// </summary>
        public MultiReadMutationScorer Clone()
        {
            // Invoke copy constructor
            return new MultiReadMutationScorer(this);
        }


        public float BaselineScore()
        {
            return scorer.BaselineScore();
        }


        /// <summary>
        /// Return the ratio of the P(T'|O) / P(T|O) where T is the template that the MutationEvaluator was
        /// constructed with, and T' is the template resulting from the application of Mutation m to template T.
        /// </summary>
        /// <param name="m">The mutation to apply to the original template</param>
        /// <returns>The ratio of the mutated and original template likelihoods</returns>
        public float ScoreMutation(Mutation m)
        {
            return scorer.Score(m.Type, m.TemplatePosition, m.Base);
        }
        
        public float ScoreMutationFast(Mutation m)
        {
            using (var ccMutation = m.ToConsensusCoreType())
            {
                return scorer.FastScore(ccMutation);
            }
        }

        public float ScoreMutationWeighted(Mutation m, float[] mappingRatios)
        {
            var scores = GetScores(m);

            var acc = 0.0f;

            for (int i = 0; i < scores.Length; i++)
            {
                if (mappingRatios[i] > 1e7)
                {
                    // Read maps essentially uniquely to this haplotype.
                    acc += scores[i];
                }
                else if (mappingRatios[i] < 1e-7)
                {
                    // Doesn't map strongy enough to affect outcome.
                }
                else
                {
                    // Partial mapping
                    var expScore = (Math.Exp(scores[i]) + (1.0f/mappingRatios[i])) / (1.0f + (1.0f/mappingRatios[i]));
                    var score = Math.Log(expScore);
                    acc += (float) score;
                }
            }

            return acc;
        }

        /// <summary>
        /// Score a compound mutation by packing the single base mutations into a single multi-base mutation.
        /// FIXME - not debugged in ConsensusCore yet.
        /// </summary>
        public float ScoreCompoundMutation(List<Mutation> mutations)
        {
            var ccMuts = mutations.Map(m => new ConsensusCore.Mutation(m.Type, m.TemplatePosition, m.Base));

            try
            {
                var mutVect = new MutationVector();
                foreach (var m in ccMuts)
                {
                    mutVect.Add(m);
                }

                return 0.0f;

                // FIXME support 'compound' mutations
                // var combinedMut = ConsensusCore.Mutation.Combine(mutVect, scorer.Template());
                // return scorer.Score(combinedMut);
            }
            finally
            {
                // Clean up C++ mutation objects
                ccMuts.ForEach(m => m.Dispose());
            }
        }

        /// <summary>
        /// Get the Quiver score deltas for this mutation for each read with a default value for unobserved mutations
        /// </summary>
        public float[] GetScoresWithDefault(Mutation m, float defaultValue)
        {
            using (var scores = scorer.Scores(m.Type, m.TemplatePosition, m.Base, defaultValue))
            {
                var r = new float[scores.Count];
                scores.CopyTo(r);
                return r;
            }
        }

        /// <summary>
        /// Aggregate the quiver score deltas according to some supplied function
        /// </summary>
        public float AggregateScoresWithDefault(Mutation m, float defaultValue, Func<float, float, float> f)
        {
            using (var scores = scorer.Scores(m.Type, m.TemplatePosition, m.Base, defaultValue))
            {
                return scores.Aggregate(defaultValue, f);
            }
        }

        /// <summary>
        /// Get the Quiver score deltas for this mutation for each read.
        /// </summary>
        public float[] GetScores(Mutation m)
        {
            return GetScoresWithDefault(m, 0);
        }

        /// <summary>
        /// Get the baseline Quiver score of each read.
        /// </summary>
        public float[] GetBaselineScores()
        {
            using (var scores = scorer.BaselineScores())
            {
                var r = new float[scores.Count];
                scores.CopyTo(r);
                return r;
            }
        }

        /// <summary>
        /// Apply a set of mutations to the template
        /// </summary>
        public void ApplyMutations(List<Mutation> mutations)
        {
            // Convert to the CC mutation type
            // NOTE - you must ensure that the CC mutations are not GC'd/deleted before the C++ ApplyMutations method is finished.
            var ccMuts =
                mutations.Map(m => new ConsensusCore.Mutation(m.Type, m.TemplatePosition, m.Base));
            try
            {
                ApplyToScorer(ccMuts);
            }
            catch (Exception e)
            {
                Console.WriteLine(e.Message);
                if(e.InnerException != null)
                    Console.WriteLine(e.InnerException.Message);
                Console.WriteLine("Got bad mutations!");
                //System.Diagnostics.Debugger.Break();
            }

            // Clean up C++ mutation objects
            ccMuts.ForEach(m => m.Dispose());
        }

        private void ApplyToScorer(IEnumerable<ConsensusCore.Mutation> muts)
        {
            var mutVect = new MutationVector();
            foreach(var m in muts)
                mutVect.Add(m);

            scorer.ApplyMutations(mutVect);
        }

        /// <summary>
        /// Current template
        /// </summary>
        public TrialTemplate Template
        {
            get { return new TrialTemplate(scorer.Template(), StartAdapterBases, EndAdapterBases); }
        }


        /// <summary>
        /// Performs application-defined tasks associated with freeing, releasing, or resetting unmanaged resources.
        /// </summary>
        public void Dispose()
        {
            if (scorer != null)
            {
                scorer.Dispose();
                scorer = null;
            }
        }

        /// <summary>
        /// Static helper method for writing out a Quiver score matrix.  Useful for making test datasets for phasing.
        /// </summary>
        public static void SaveMultiSplitScores(MultiReadMutationScorer sc, string filename)
        {
            // Use all possible substitutions as potential splitting mutations
            var muts =
                GenerateMutations.GenerateUniqueMutations(sc.Template)
                    .Where(s => s.Type == MutationType.SUBSTITUTION)
                    .ToArray();

            var mutScoreVectors = muts.Map(sc.GetScores);
            var positions = muts.Map(m => m.TemplatePosition);

            var binaryFormatter = new BinaryFormatter();
            using (
                var serializationStream = new FileStream(filename, FileMode.Create,
                    FileAccess.Write,
                    FileShare.None))
            {
                binaryFormatter.Serialize(serializationStream, new Tuple<float[][], int[]>(mutScoreVectors, positions));
            }
        }
    }
}


