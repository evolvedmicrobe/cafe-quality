#region Copyright (c) 2010, Pacific Biosciences of California, Inc.
//
// All rights reserved.
//
// THIS SOFTWARE CONSTITUTES AND EMBODIES PACIFIC BIOSCIENCES’ CONFIDENTIAL
// AND PROPRIETARY INFORMATION.
//
// Disclosure, redistribution and use of this software is subject to the
// terms and conditions of the applicable written agreement(s) between you
// and Pacific Biosciences, where “you” refers to you or your company or
// organization, as applicable.  Any other disclosure, redistribution or
// use is prohibited.
//
// THIS SOFTWARE IS PROVIDED BY PACIFIC BIOSCIENCES AND ITS CONTRIBUTORS "AS
// IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
#endregion

using System;
using System.Collections.Generic;
using System.Linq;
using ConsensusCore;
using PacBio.IO;
using PacBio.Utils;
using PacBio.Align;

namespace PacBio.Consensus
{
    /// <summary>
    /// Represents the current guess at the template for a single molecule consensus procedure 
    /// </summary>
    public class TrialTemplate
    {
        public TrialTemplate()
        {
            
        }

        public TrialTemplate(string template, int startAdapterBases = 0, int endAdapterBases = 0)
        {
            Sequence = template;
            StartAdapterBases = startAdapterBases;
            EndAdapterBases = endAdapterBases;
        }

        /// <summary>
        /// Current guess at the template sequence, expressed in the forward strand
        /// </summary>
        public string Sequence;

        /// <summary>
        /// Number of bases of the adapter included at front of sequence
        /// </summary>
        public int StartAdapterBases;

        /// <summary>
        /// Number of bases of adapter included at the end of sequence
        /// </summary>
        public int EndAdapterBases;

        /// <summary>
        /// Length of the template
        /// </summary>
        public int Length
        {
            get { return Sequence.Length; }
        }

        /// <summary>
        /// Get the template sequence in channel space for the current TrialTemplate, converted to the passed strand.
        /// </summary>
        /// <param name="st">The strand to convert the template sequence to.</param>
        /// <returns>A channel-space (0-3) template sequence</returns>
        public string GetSequence(Strand st)
        {
            switch (st)
            {
                case Strand.Forward:
                    return Sequence;

                case Strand.Reverse:
                    return Sequence.ReverseComplement();

                default:
                    throw new ArgumentException("Unrecognized Strand type");
            }
        }

        /// <summary>
        /// Set the template sequence
        /// </summary>
        /// <param name="sequence">The channel-space template sequence</param>
        /// <param name="strand">The strand that the sequence is provided in</param>
        public void SetSequence(string sequence, Strand strand)
        {
            if (strand == Strand.Forward)
            {
                Sequence = sequence;
            }
            else
            {
                Sequence = DNA.ReverseComplement(sequence);
            }
        }

        /// <summary>
        /// The channel-space sequence of the forward strand, excluding the flanking adapter regions
        /// </summary>
        public string InsertSequence()
        {
            return Sequence.Substring(StartAdapterBases, Sequence.Length - StartAdapterBases - EndAdapterBases);
        }

        /// <summary>
        /// Return a new TrialTemplate, with the template mutation m applied
        /// </summary>
        /// <param name="m">The mutation to apply to the template</param>
        /// <returns>A fresh TrialTemplate struct</returns>
        public TrialTemplate Mutate(Mutation m)
        {
            return new TrialTemplate
                       {
                           EndAdapterBases = EndAdapterBases,
                           Sequence = m.Apply(Sequence),
                           StartAdapterBases = StartAdapterBases,
                       };
        }

        /// <summary>
        /// Simultaneously apply a set of mutations to the TrialTemplate and return a fresh TrialTemplate
        /// </summary>
        /// <param name="muts">A sequence of mutations</param>
        /// <returns>A fresh TrialTemplate with the mutations applied</returns>
        public TrialTemplate Mutate(IEnumerable<Mutation> muts)
        {
            return new TrialTemplate
            {
                EndAdapterBases = EndAdapterBases,
                Sequence = Mutation.ApplyMany(muts.ToList(), Sequence),
                StartAdapterBases = StartAdapterBases
            };
        }
    }

    /// <summary>
    /// A simple struct for tracking how beneficial a mutation is to the CRF likelihood.  
    /// When creating a MutationScore, set it's Exists field to true, so that we can use struct semantics and have an existence check.
    /// </summary>
    public struct MutationScore
    {
        public Mutation Mutation;
        public double Score;
        public bool Exists;
    }

    /// <summary>
    /// A collection of methods for finding the consensus of bunch of reads from their pulse data using
    /// a probabilistic alignment model.
    /// </summary>
    public class FindConsensus
    {
        /// <summary>
        /// CCS POA implementation. NOTE! This method is a bit hairy, so tread with caution.
        /// The basic tasks are:
        ///     1. Pick an ordering of the subreads to feed to POA
        ///     2. Feed the subreads to the POA
        ///     3. Tack on the last adapterPadBases of the adapter to the template
        ///     4. Find the right read / template intervals, including the adapter padding
        ///     5. Also implement 'sticky ends' to push the read alignment to the end of the template where possible.
        /// </summary>
        public static TrialTemplate InitialConsensusTemplate(
            RegionSets regionSets, IZmwBases bases, out float consensusScore,
            out List<AlignedSequenceReg> regions, int adapterPadBases, string adapter)
        {
            //var insertRegions = regionSets.InsertRegions;
            var sequence = bases.Sequence;
            var qvs = bases.QV;

            var avgLength = regionSets.InsertRegions.Select(r => r.Length).Average();
            Func<int, double> lengthPenalty = l => Math.Min((double)l / avgLength, (double)avgLength / l);

            // convert subreads into strings on forward strand.
            var subReadOrder = regionSets.InsertRegions.Select((reg, idx) => 
            {
                var l = Math.Max(0, reg.Length);
                var rQv = qvs.Slice(reg.Start, l);

                var predAcc = 0.0;        

                if (reg.Length > 0)
                    predAcc = (1.0 - rQv.Select(QVs.QVToProb).Average()) * lengthPenalty(reg.Length);

                return new { Acc = predAcc, Idx = idx };
            }).OrderByDescending(v => v.Acc).Select(v => v.Idx).ToArray();

            var subreads = regionSets.InsertRegions.Select(reg =>
            {
                var l = Math.Max(0, reg.Length);
                string insertSeqBases = sequence.Substring(reg.Start, l);

                if (reg.Strand == Strand.Forward)
                    return insertSeqBases;
                else
                    return insertSeqBases.ReverseComplement();
            }).ToArray();

            subreads = subReadOrder.Map(v => subreads[v]);
            var insertRegions = subReadOrder.Map(v => regionSets.InsertRegions[v]);
            var paddedRegions = subReadOrder.Map(v => regionSets.PaddedInsertRegions[v]);

            // put the reads into the 'local' poa
            var b = new PoaLocal();
            b.AddReads(subreads);

            // Get out the POA consensus and the extents of each reads alignment to the consensus
            string initialTemplateSequence;
            var readExtents = b.FindConsensusAndAlignments(1, out consensusScore, out initialTemplateSequence);


            // Combine the POA extents and the raw read extents of each subread into the final extents to use in consensus
            //  - Use 'sticky-ends' -- if the POA extent comes within n bases of an adapter hit, just included the adapter hit.
            //  - Use 'roll-in'  for ends without an adapter hit, we will use 'unpinned' boundary conditions in the recursion, so the 
            //    template bounds should be pushed back k bases from the observed extent.
            //  - Keep track of which ends actually have an adapter hit
            var stickyLength = 15;

            var padAtSides = insertRegions.Select((reg, i) =>
            {
                // Some subreads may not hit the POA at all -- remove them.
                if (!readExtents.ContainsKey(i))
                    return null;

                var extents = readExtents[i];

                // Don't bother with short subreads
                if (extents.Item2.Read - extents.Item1.Read < 15)
                    return null;

                var hitLeft = false;
                var hitRight = false;

                var adaLeft = reg.Strand == Strand.Forward ? reg.AdapterHitBefore : reg.AdapterHitAfter;
                var adaRight = reg.Strand == Strand.Forward ? reg.AdapterHitAfter : reg.AdapterHitBefore;

                // Sticky check front
                if (adaLeft && extents.Item1.Template < stickyLength && extents.Item1.Read < stickyLength)
                {
                    hitLeft = true;
                }

                // Sticky check back
                if (adaRight && (initialTemplateSequence.Length - extents.Item2.Template) < stickyLength && (reg.Length - extents.Item2.Read) < stickyLength)
                {
                    hitRight = true;
                }

				return new Tuple<bool, bool>(hitLeft, hitRight);

            }).Where(v => v != null).ToList();

            var adaPadLeft = padAtSides.Any(v => v.Item1);
            var adaPadRight = padAtSides.Any(v => v.Item2);

            var startOffset = adaPadLeft ? adapterPadBases : 0;
            
            var ttpl = new TrialTemplate()
            {
                EndAdapterBases = 0,
                StartAdapterBases = 0
            };

            // Tack the adapters onto the POA sequence if any of the reads get to that end of the molecule:
            // i.e. Only pad on ends where there is at least 1 adapter hit.
            var fullSeq = initialTemplateSequence;

            if (adaPadLeft)
            {
                var adaLength = adapter.Length;
                var adapterEnd = adapter.Substring(adaLength - adapterPadBases, adapterPadBases);
                fullSeq = adapterEnd + fullSeq;

                ttpl.StartAdapterBases = adapterPadBases;
            }

            if (adaPadRight)
            {
                var adapterStart = adapter.Substring(0, adapterPadBases);
                fullSeq = fullSeq + adapterStart;

                ttpl.EndAdapterBases = adapterPadBases;
            }

            ttpl.SetSequence(fullSeq, Strand.Forward);


            regions = insertRegions.Select((reg, i) =>
            {
                // Some subreads may not hit the POA at all -- remove them.
                if (!readExtents.ContainsKey(i))
                    return null;

                var extents = readExtents[i];

                // Don't bother with short subreads
                if (extents.Item2.Read - extents.Item1.Read < 15)
                    return null;

                if (reg.Strand == Strand.Forward)
                {
                    var rFront = reg.Start + extents.Item1.Read;
                    var tFront = extents.Item1.Template + startOffset;

                    var rEnd = reg.Start + extents.Item2.Read;
                    var tEnd = extents.Item2.Template + startOffset;

                    var hitBefore = false;
                    var hitAfter = false;

                    // Sticky check front
                    if (reg.AdapterHitBefore && extents.Item1.Template < stickyLength && extents.Item1.Read < stickyLength)
                    {
                        var paddedRegion = paddedRegions[i];
                        rFront = paddedRegion.Start;
                        tFront = 0;
                        hitBefore = true;
                    }

                    // Sticky check back
                    if (reg.AdapterHitAfter && (initialTemplateSequence.Length - extents.Item2.Template) < stickyLength && (reg.Length - extents.Item2.Read) < stickyLength)
                    {
                        var paddedRegion = paddedRegions[i];
                        rEnd = paddedRegion.End;
                        tEnd = fullSeq.Length;
                        hitAfter = true;
                    }

                    var alReg = new AlignedSequenceReg(rFront, rEnd, tFront, tEnd, reg.Strand)
                    {
                        AdapterHitBefore = hitBefore,
                        AdapterHitAfter = hitAfter
                    };
                    return alReg;
                }
                else
                {
                    var rFront = reg.Start + (reg.Length - extents.Item2.Read - 1);
                    var tFront = extents.Item1.Template + startOffset;

                    var rEnd = reg.Start + (reg.Length - extents.Item1.Read - 1);
                    var tEnd = extents.Item2.Template + startOffset;

                    var hitBefore = false;
                    var hitAfter = false;

                    // Sticky check front
                    //if (reg.AdapterHitAfter && extents.Item1.Template < stickyLength)
                    if (reg.AdapterHitAfter && extents.Item1.Template < stickyLength && extents.Item1.Read < stickyLength)
                    {
                        var paddedRegion = paddedRegions[i];
                        rEnd = paddedRegion.End;
                        tFront = 0;
                        hitAfter = true;
                    }

                    // Sticky check back
                    //if (reg.AdapterHitBefore && (initialTemplateSequence.Length - extents.Item2.Template) < stickyLength)
                    if (reg.AdapterHitBefore && (initialTemplateSequence.Length - extents.Item2.Template) < stickyLength && (reg.Length - extents.Item2.Read) < stickyLength)
                    {
                        var paddedRegion = paddedRegions[i];
                        rFront = paddedRegion.Start;
                        tEnd = fullSeq.Length;
                        hitBefore = true;
                    }

                    var alReg = new AlignedSequenceReg(rFront, rEnd, tFront, tEnd, reg.Strand)
                    {
                        AdapterHitBefore = hitBefore,
                        AdapterHitAfter = hitAfter
                    };
                    return alReg;
                }
            }).Where(v => v != null).ToList();


            return ttpl;
        }

        /// <summary>
        /// Main Quiver loop and QV estimation
        /// </summary>
        public static IZmwConsensusBases MultiReadConsensusAndQv(
            MultiReadMutationScorer scorer, List<AlignedSequenceReg> regions, 
            ISequencingZmw zmw, int iterations, out int iterationsTaken)
        {
            int numberOfAcceptedMutations;
            int numberOfTriedMutations;
            // Run the quiver refinement step
            var result = InnerMultiReadConsensus(scorer, iterations, out iterationsTaken, out numberOfAcceptedMutations, out numberOfTriedMutations);

            // Sort the regions that will get reported in the bas.h5
            var sortedRegions = regions.Map(v => (DelimitedSeqReg) v).OrderBy(v => v.Start).ToArray();

            // Compute the QVs from the mutation scores
            var qvData = ComputeConsensusQs(result.Item2, result.Item1);
            var res = ComputeConsensusQVs(sortedRegions, zmw, qvData);
            // Reset so that that the number of passes is the number
            // of scored subreads, not the insert regions (can differ due to alpha/beta mismatch).
            res.NumPasses = scorer.GetBaselineScores ().Length;
            res.NumberOfMutations = numberOfAcceptedMutations;
            res.NumberOfTriedMutations = numberOfTriedMutations;

            return res;
        }

        /// <summary>
        /// Quiver refinement loop
        /// </summary>
        public static Tuple<TrialTemplate, List<MutationScore>> InnerMultiReadConsensus(MultiReadMutationScorer scorer,
                                                                                        int iterations,
                                                                                        out int iterationsTaken, 
                                                                                        out int mutationsAccepted,
                                                                                        out int mutationsTested)
        {
            mutationsAccepted = 0;
            mutationsTested = 0;
            Func<Mutation, MutationScore> scoreMutation = 
                m => new MutationScore { Score = scorer.ScoreMutationFast(m), Mutation = m, Exists = true };


            // This will be used to derive QVs
            List<MutationScore> allScores = null;
            int phase = 0;
            double score;

            int mutationSpacing = 9;

            // This corresponds to a probability of ~51% - to avoid flip-flopping between opposite mutations
            float minScore = 0.35f;
            int prevMutationWindow = 12;

            Func<IEnumerable<Mutation>, List<Mutation>> screenMutations = mutationsToTry => FindMutations(mutationsToTry, scoreMutation, out score, mutationSpacing, minScore);

            var tpl = scorer.Template;

            // This will track the last batch of mutations we made
            List<Mutation> muts = null;
            var i = 0;

            Action<List<Mutation>> setMuts = m => { muts = m; };


            // there are a few phases here where we attempt different types of mutations in different regio
            while (true)
            {
                tpl = scorer.Template;

                // If we are over our iteration allotment, compute the QVs and bail
                if (i >= iterations)
                {
                    phase = 2;
                }

                var p = phase;
                switch (p)
                {

                    case 0:
                        // Try all mutations - this phase only runs once, then we'll only test mutations close to the ones we found
                        var mutationsToTry = GenerateMutations.GenerateUniqueMutations (tpl, false).ToList();
                        mutationsTested += mutationsToTry.Count;
                        setMuts(screenMutations(mutationsToTry));
                        phase = 1;

                        if (muts.Count == 0)
                        {
                            // if we don't get any mutations here, go to the end
                            phase = 2;
                        }

                        break;

                        
                case 1:
                        // Test mutations near previously applied batch, until no new ones are found.
                        // Only test for subs when we are getting toward the end.
                    var testSubs = i > 4;
                    mutationsToTry = GenerateMutations.PrevMutationFilter (GenerateMutations.GenerateUniqueMutations (tpl, testSubs), muts, prevMutationWindow).ToList ();
                    mutationsTested += mutationsToTry.Count;
                        setMuts(screenMutations(mutationsToTry));

                        if (muts.Count == 0)
                            phase = 2;

                        break;

                    case 2:
                        // In phase 3 we shouldn't get many new mutations. Compute the scores
                        // of all possible mutations in preparation for computing the QVs.
                        // If we do find real mutations, apply them and try again.
                        allScores = UniqueMutationsScores(scoreMutation, tpl);
                        var possibleMutations =
                            allScores.Where(s => !s.Mutation.IsSynonymous(tpl) && s.Score > minScore).ToList();
                        var bestMutations = SpacedSelector.BestMutations(possibleMutations, mutationSpacing);
                        setMuts(bestMutations.Select(m => m.Mutation).ToList());
                        break;
                }

                // If we get to the max iterations, or have a small number of mutations, then bail
                if ((allScores != null) && (i >= iterations || muts.Count == 0))
                {
                    break;
                }

                // Apply the current set of mutations. If there are none, then we
                // will save some time by avoid the creation of new F/B matrices.
                if (muts.Count > 0) {
                    mutationsAccepted += muts.Count;
                    scorer.ApplyMutations (muts);
                }
                i++;
            }
            //var lastScore = scorer.BaselineScore ();
            /* var mutations = new List<Mutation> { new Mutation () {
                    TemplatePosition = 42,
                    Type = MutationType.DELETION,
                    Base = '-'
                },
                new Mutation () { TemplatePosition = 38, Type = MutationType.INSERTION, Base = 'G' }
            };
            scorer.ApplyMutations (mutations);
            var newScore = scorer.BaselineScore ();
            Console.WriteLine (lastScore - newScore);
*/
            iterationsTaken = i;
            return new Tuple<TrialTemplate, List<MutationScore>>(tpl, allScores);
        }


        public static int ImproveConsensus(MultiReadMutationScorer scorer, int maxIterations, double[] mappingRatios = null, bool fast = false)
        {
            Func<Mutation, MutationScore> scoreMutation;

            if (mappingRatios == null)
            {
                if (fast)
                {
                    scoreMutation = m => new MutationScore { Score = scorer.ScoreMutationFast(m), Mutation = m, Exists = true };
                }
                else
                {
                    // Normal style -- only one template in play
                    scoreMutation = m => new MutationScore { Score = scorer.ScoreMutation(m), Mutation = m, Exists = true };
                }
            }
            else
            {
                // Multi-template scoring -- the mapping ratio is P_this_read / sum_i (P_read_i)
                scoreMutation = m => new MutationScore { Score = scorer.ScoreMutationWeighted(m, mappingRatios), Mutation = m, Exists = true };
            }

            double score;

            // 0.35f corresponds to a probability of ~51% - to avoid flip-flopping between opposite mutations
            // lhepler: I don't know how Pat computed that number, I get 41.3%
            const float minScore = 0.35f;
            const int mutationSpacing = 8;
            const int prevMutationWindow = 22;

            Func<IEnumerable<Mutation>, List<Mutation>> screenMutations =
                mutationsToTry => FindMutations(mutationsToTry, scoreMutation, out score, mutationSpacing, minScore);
                    // FIXME -- consider the compound mutation search procedure
                    //FindMutationsAndSearch(mutationsToTry, scoreMutation, mutationSpacing, minScore, compoundScore);

            var tplHistory = new HashSet<string>();
            var tpl = scorer.Template;
            int iters;

            tplHistory.Add(tpl.GetSequence(Strand.Forward));

            // This will track the last batch of mutations we made, try all mutations
            var muts = screenMutations(GenerateMutations.GenerateUniqueMutations(tpl, true));

            for (iters = 0; muts.Any() && iters < maxIterations; iters++)
            {
                // Apply the current set of mutations. If there are none, then we
                // will save some time by avoid the creation of new F/B matrices.

                // cycle avoidance
                if (muts.Count > 1)
                {
                    var nextTpl = tpl.Mutate(muts);

                    if (tplHistory.Contains(nextTpl.GetSequence(Strand.Forward)))
                    {
                        muts = muts.Take(1).ToList();
                    }
                }

                scorer.ApplyMutations(muts);

                tpl = scorer.Template;

                tplHistory.Add(tpl.GetSequence(Strand.Forward));

                var mutationsToTry = GenerateMutations.PrevMutationFilter(GenerateMutations.GenerateUniqueMutations(tpl, true), muts, prevMutationWindow);
                muts = screenMutations(mutationsToTry);
            }

            return iters;
        }

        public static List<MutationScore> ComputeAllQVs(MultiReadMutationScorer scorer, double[] mappingRatios = null)
        {
            Func<Mutation, MutationScore> scoreMutation;

            if (mappingRatios == null)
            {  
                // Normal style -- only one template in play
                scoreMutation = m => new MutationScore { Score = scorer.ScoreMutation(m), Mutation = m, Exists = true };
            }
            else
            {
                // Multi-template scoring -- the mapping ratio is P_this_read / sum_i (P_read_i)
                scoreMutation = m => new MutationScore { Score = scorer.ScoreMutationWeighted(m, mappingRatios), Mutation = m, Exists = true };
            }
            
            var tpl = scorer.Template;
            var allScores = UniqueMutationsScores(scoreMutation, tpl);
            return allScores;
        }

        private static double MutationProbability(MutationScore sc)
        {
            if (!sc.Exists)
                return 0.0;

            return ScoreToErrorProb(sc.Score);
        }

        /// <summary>
        /// Mutation scores are logit transformed measured of the probability the that mutation is correct
        /// </summary>
        /// <param name="score"></param>
        /// <returns></returns>
        public static double ScoreToErrorProb(double score)
        {
            // Use a logit transform
            var errProb = (1.0 / (1.0 + Math.Exp(-score)));

            // Cap the apparent QV at 50
            return errProb;
        }

        public class ConsensusQVs
        {
            /// <summary>
            /// Predicted accuracy of called sequence
            /// </summary>
            public float PredictedAccuracy { get; set; }

            /// <summary>
            /// Phred-like quality value of the base
            /// </summary>
            public IList<byte> QV { get; set; }

            /// <summary>
            /// The called base sequence
            /// </summary>
            public string Sequence { get; set; }

            /// <summary>
            /// Phred-style quality value indicating the probability that the current base is an insertion
            /// </summary>
            public IList<byte> InsertionQV { get; set; }

            /// <summary>
            /// Phred-style quality value indicating the total probability of a deleted base before the current base,
            /// excluding the merge case.
            /// </summary>
            public IList<byte> DeletionQV { get; set; }

            /// <summary>
            /// Phred-style quality value indicating the total probability that a merged pulse-call produced this base
            /// </summary>
            public IList<byte> MergeQV { get; set; }

            /// <summary>
            /// The ASCII code of the most likely base to have been deleted before the current base.
            /// The ASCII code for 'N' represents a 'roughly constant' deletion probability over all bases. 
            /// </summary>
            public IList<char> DeletionTag { get; set; }

            /// <summary>
            /// Phred-style quality value indicating the total probability that the current basecall is a substitution error
            /// </summary>
            public IList<byte> SubstitutionQV { get; set; }

            /// <summary>
            /// The ASCII code of the most likely alternative basecall at this position
            /// </summary>
            public IList<char> SubstitutionTag { get; set; }
        }


        /// <summary>
        /// Compute consensus QVs by summarizing the mutation scores for each position
        /// </summary>
        public static ZmwConsensusBases ComputeConsensusQVs(IList<DelimitedSeqReg> regions,  ISequencingZmw zmw, ConsensusQVs qvData)
        {

            var zmwBases = new ZmwConsensusBases(zmw, regions, qvData.PredictedAccuracy, qvData.Sequence.Length, -999)
                               {
                                   Base = qvData.Sequence.ToCharArray(),
                                   InsertionQV = qvData.InsertionQV,
                                   DeletionQV = qvData.DeletionQV,
                                   DeletionTag = qvData.DeletionTag,
                                   PulseIndex = qvData.Sequence.Length.Fill(i => 0),
                                   QV = qvData.QV,
                                   SubstitutionQV = qvData.SubstitutionQV,
                                   SubstitutionTag = qvData.SubstitutionTag,
                               };
            return zmwBases;
        }

        public static ConsensusQVs ComputeConsensusQs(List<MutationScore> score, TrialTemplate tpl)
        {
            var pad = tpl.StartAdapterBases;
            var insert = tpl.InsertSequence();

            Func<MutationScore, int> idx = m => m.Mutation.TemplatePosition - pad;

            var n = insert.Length;
            var deleteScores = new MutationScore[insert.Length];
            var insertScores = new MutationScore[insert.Length,4];
            var substitutionScores = new MutationScore[insert.Length,4];

            //var icare = score.Where (x => x.Mutation.TemplatePosition > 37 && x.Mutation.TemplatePosition < 42).ToList();
            //score.Sort ((x, y) => -x.Score.CompareTo (y.Score));
            //Console.WriteLine (icare);

            // Distribute the Mutation scores by mutation type and base, aligned with the insert sequence
            foreach (var msc in score)
            {
                var i = idx(msc);

                if (i < n)
                {
                    switch (msc.Mutation.Type)
                    {
                        case MutationType.DELETION:
                            deleteScores[i] = msc;
                            break;

                        case MutationType.INSERTION:
                            insertScores[i, DNA.BaseCode[msc.Mutation.Base]] = msc;
                            break;

                        case MutationType.SUBSTITUTION:
                            if (msc.Mutation.Base != insert[msc.Mutation.TemplatePosition - pad])
                            {
                                substitutionScores[i, DNA.BaseCode[msc.Mutation.Base]] = msc;
                            }
                            break;
                    }
                }
            }

            // Compute insertion 
            var insertQv = deleteScores.Map(v => QVs.ProbToQV(MutationProbability(v)));
            
            // DeletionQV for each base at each position
            var perBaseDelProb = insertScores.Map(MutationProbability);
            var delTag = perBaseDelProb.Collapse(scs => DNA.Bases[scs.IMax()], 1);

            // Sum up the deletion QVs at each position
			var totalDelQv = perBaseDelProb.Collapse(scs => QVs.ProbToQV(QVs.CombineErrorProbability(scs)), 1);

            // Substitution QV for each base at each position
            var perBaseSubsProb = substitutionScores.Map(MutationProbability);
            var subsTag = perBaseSubsProb.Collapse(scs => DNA.Bases[scs.IMax()], 1);

            // Sum up substitution QV at each position
			var totalSubsQv = perBaseSubsProb.Collapse(scs => QVs.ProbToQV(QVs.CombineErrorProbability(scs)), 1);

            var totalQv = insert.Length.Fill(
                i =>
                {
                    var v = QVs.CombineErrorProbability(
                        MutationProbability(deleteScores[i]),
                        perBaseDelProb[i, 0], perBaseDelProb[i, 1],
                        perBaseDelProb[i, 2], perBaseDelProb[i, 3],
                        perBaseSubsProb[i, 0],perBaseSubsProb[i, 1],
                        perBaseSubsProb[i, 2],perBaseSubsProb[i, 3]);

                    return QVs.ProbToQV(v);
                });

            return new ConsensusQVs
                {
                    Sequence = insert,
                    InsertionQV = insertQv,
                    DeletionQV = totalDelQv,
                    DeletionTag = delTag,
                    QV = totalQv,
                    SubstitutionQV = totalSubsQv,
                    SubstitutionTag = subsTag,
                    PredictedAccuracy =  1.0f - (float) totalQv.Select(QVs.QVToProb).Average()
                };
        }
        
        /// <summary>
        /// Generate a list of beneficial template mutations, given an initial trial template, a series of PulsePassModels, and some
        /// constraints on the density and score of the returned mutations
        /// </summary>
        public static List<Mutation> FindMutations(IEnumerable<Mutation> mutations, 
            Func<Mutation, MutationScore> scoreMutation, out double score, int minSpacing, float minScore)
        {
            // Generate all possible mutations of the template.  GenerateAllMutations will not generate mutations
            // in the adapter regions. Filter out mutations that don't meet the minimum score threshold.
            var possibleMutations = mutations.Select(scoreMutation).Where(ms => ms.Score > minScore).ToList();

            // Find the set of mutations with at minSpacing template bases between mutations with the highest total score
            var bestMutations = SpacedSelector.BestMutations(possibleMutations, minSpacing);
            
            // Get the total score of the selected mutations
            score = bestMutations.Select(ms => ms.Score).Sum();
            // Return the muatations
            var selectedMutations = bestMutations.Select(s => s.Mutation).ToList();
            return selectedMutations;
        }

        /// <summary>
        /// Generate a list of beneficial template mutations, given an initial trial template, a series of PulsePassModels, and some
        /// constraints on the density and score of the returned mutations
        /// </summary>
        static List<MutationScore> UniqueMutationsScores(Func<Mutation, MutationScore> scoreMutation, TrialTemplate tpl)
        {
            // Generate all possible mutations of the template.  GenerateAllMutations will not generate mutations
            // in the adapter regions. Filter out mutations that don't meet the minimum score threshold.
            return GenerateMutations.GenerateUniqueMutations(tpl).Select(
                m =>
                    {
                        if (m.IsSynonymous(tpl))
                            return new MutationScore
                                {
                                    Exists = true,
                                    Mutation = m,
                                    Score = 0
                                };
                        else
                            return scoreMutation(m);
                    }).ToList();
        }
    }
}