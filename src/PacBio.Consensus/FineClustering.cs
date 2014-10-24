using System;
using System.Collections.Generic;
using System.Linq;
using ConsensusCore;
using PacBio.Align;
using PacBio.IO;
using PacBio.Utils;
using PacBio.FSharp.Utils;

namespace PacBio.Consensus
{
    /// <summary>
    /// 'Fine clustering' is a synonym for recursive phasing. We get a bunch of subread that are known to come from a collection of similar template
    /// sequences, and we want to estimate the template sequences and which reads came from which templates.
    /// </summary>
    public class FineClustering
    {
        public static bool DoPhasing = true;

        private static PacBioLogger logger = PacBioLogger.GetLogger("FineClustering");

        private static void Log(LogLevel level, string msg)
        {
            logger.Log(level, msg);
        }

        private static void Log(LogLevel level, string msg, params object[] args)
        {
            logger.Log(level, String.Format(msg, args));
        }


        public static Tuple<string, AlignedSequenceReg, IZmwBases> GetMappedExtents(Subread r, string template)
        {
            var seq = r.FwdSequence;
            var rseq = seq.ReverseComplement();

            // Do a sparse alignment of the reads to get the alignment band
            var sparse = new SparseAligner(6);

            var fwdFragments = sparse.SparseAlignStrings(seq, template);
            var revFragments = sparse.SparseAlignStrings(rseq, template);

            Strand strand;
            List<Fragment> fragments;

            if (fwdFragments.Count > revFragments.Count)
            {
                strand = Strand.Forward;
                fragments = fwdFragments;
            }
            else
            {
                strand = Strand.Reverse;
                fragments = revFragments;
                seq = rseq;
            }

            SparseBand band;
            
            try
            {
                band = new SparseBand(fragments, 12, seq.Length, template.Length);
            }
            catch (SparseAlignmentException)
            {
                return null;
            }

            var cells = BandedAlignment.Local(seq, template, band);

            var subReadStart = r.Region.Start;
            var subReadLength = r.Region.Length;

            var startPt = cells.First();
            var endPt = cells.Last();

            AlignedSequenceReg res;

            if (strand == Strand.Forward)
            {
                res = new AlignedSequenceReg(
                    subReadStart + startPt.Read,
                    subReadStart + (endPt.Read + 1),
                    startPt.Template,
                    endPt.Template + 1,
                    strand);
            }
            else
            {
                res = new AlignedSequenceReg(
                    subReadStart + subReadLength - (endPt.Read + 1),
                    subReadStart + subReadLength - startPt.Read,
                    startPt.Template,
                    endPt.Template + 1,
                    strand);
            }

            res.Accuracy = (float)AlignCell.Accuracy(cells);

            return new Tuple<string, AlignedSequenceReg, IZmwBases>(r.SubreadId, res, r.Bases);
        }


        public static string GetConsensusAndOverlaps(Subread[] reads, out Tuple<AlignedSequenceReg, IZmwBases>[] extents)
        {
            const int cleanThreshold = 12;

            // Stuff them into a POA
            var poa = new PoaLocal();

            var matchedStrand = new Dictionary<int, Strand>();
            var addedSubread = new Dictionary<int, Subread>();
            int i = 0;

            // Make the POA
            foreach (var r in reads)
            {
                var seq = r.FwdSequence;
                var rseq = seq.ReverseComplement();

                if (poa.NumReads == 0)
                {
                    poa.AddReadSparse2(seq);
                    matchedStrand[i] = Strand.Forward;
                    addedSubread[i] = r;
                    i++;
                    continue;
                }

                using (var sFwd = poa.TryAddRead(seq))
                using (var sRev = poa.TryAddRead(rseq))
                {
                    // 70% of the current sequence or the avg length of the sequences already in the POA
                    var minScoreToAddRead = 0.7 * Math.Min(seq.Length, poa.ReadLengths.Average());

                    if (sFwd.Score > sRev.Score && sFwd.Score > minScoreToAddRead)
                    {
                        sFwd.Commit();
                        matchedStrand[i] = Strand.Forward;
                        addedSubread[i] = r;
                        i++;
                    }
                    else if (sRev.Score > sFwd.Score && sRev.Score > minScoreToAddRead)
                    {
                        sRev.Commit();
                        matchedStrand[i] = Strand.Reverse;
                        addedSubread[i] = r;
                        i++;
                    }
                }

                // Clean out crappy POA nodes
                if (poa.NumReads > cleanThreshold && poa.NumReads%cleanThreshold == 0)
                {
                    poa.ComputeCoverage();
                    poa.Graph.RemoveVertexIf(v => v.Coverage == 1);
                }
            }

            float score;
            string consensus;
            var alignments = poa.FindConsensusAndAlignments(1, out score, out consensus);
            
            extents = alignments.Select(a =>
            {
                var subReadId = a.Key;
                
                var subRead = addedSubread[subReadId];
                var strand = matchedStrand[subReadId];

                var subReadStart = subRead.Region.Start;
                var subReadLength = subRead.Region.Length;

                var startPt = a.Value.Item1;
                var endPt = a.Value.Item2;

                AlignedSequenceReg r;

                if (matchedStrand[subReadId] == Strand.Forward)
                {
                    r = new AlignedSequenceReg(
                        subReadStart + startPt.Read,
                        subReadStart + (endPt.Read + 1),
                        startPt.Template,
                        endPt.Template + 1,
                        strand);
                }
                else
                {
                    r = new AlignedSequenceReg(
                        subReadStart + subReadLength - (endPt.Read + 1),
                        subReadStart + subReadLength - startPt.Read,
                        startPt.Template,
                        endPt.Template + 1,
                        strand);
                }
                
                return new Tuple<AlignedSequenceReg, IZmwBases>(r, subRead.Bases);
            }).ToArray();

            return consensus;
        }


        private static ConsensusResult[] LabelDuplicates(ConsensusResult[] results)
        {
            for (int i = 0; i < results.Length; i++)
            {
                // if one of the sequences is a duplicate already, skip
                if (!String.IsNullOrEmpty(results[i].DuplicateOf))
                    continue;

                var iSeq = results[i].Sequence;

                for (int j = i + 1; j < results.Length; j++)
                {
                    if (!String.IsNullOrEmpty(results[j].DuplicateOf))
                        continue;

                    var jSeq = results[j].Sequence;

                    // if j is a subset of i, then it's a duplicate
                    // OR if i is a subset of j, then it's a duplicate and done
                    if (iSeq.IndexOf(jSeq) >= 0 || iSeq.IndexOf(jSeq.ReverseComplement()) >= 0)
                    {
                        results[j].DuplicateOf = results[i].FastaName;
                    }
                    else if (jSeq.IndexOf(iSeq) >= 0 || jSeq.IndexOf(iSeq.ReverseComplement()) >= 0)
                    {
                        results[i].DuplicateOf = results[j].FastaName;
                    }
                }
            }

            return results;
        }


        public static int PoaReads = 30;

        /// <summary>
        /// Main entry point for phasing and consensus.  All the subreads passed in here should map reasonably well to a common POA sequence.
        /// </summary>
        /// <param name="coarseClusterNumber">Tracked forward for logging/output purposes</param>
        /// <param name="reads">Subreads to phas</param>
        /// <param name="config">ScorerConfig of model parameters and settings</param>
        /// <param name="barcodeName">Tracked forward for logging/output purposes</param>
        /// <param name="phasingReads">Number of reads to use for phasing / consensus</param>
        /// <param name="sortByCoverage">sort the results by descending number of supporting reads</param>
        public static ConsensusResult[] FineClusteringAndConsensus(int coarseClusterNumber, Subread[] reads,
                                                                   ScorerConfig config,
                                                                   string barcodeName = "0",
                                                                   int phasingReads = 500,
                                                                   bool sortByCoverage = false)
        {
            // Compute the POA sequence using the first few reads.
            Tuple<AlignedSequenceReg, IZmwBases>[] poaExtents;
            var nPoa = Math.Min(reads.Length, PoaReads);
            var poaReads = reads.Take(nPoa).ToArray();

            // Any truncations arise from a problem in here -- POA 
            var consensus = GetConsensusAndOverlaps(poaReads, out poaExtents);

            // Map all the reads to the POA sequence to get the template extents for ConsensusCore
            var readExtents = reads
                .Select(sr => GetMappedExtents(sr, consensus))
                .Where(re => re != null)
                .ReservoirSample(phasingReads);

            if (readExtents.Length < MultiTemplateConsensus.minSplitReads)
            {
                Log(LogLevel.WARN, "Barcode '{0}', Cluster {1} -- Phasing/Consensus aborted due to low number of reads: {2} reads", barcodeName, coarseClusterNumber, readExtents.Length);
                return new ConsensusResult[] { };
            }

            Log(LogLevel.INFO, "Barcode '{0}', Cluster {1} -- Phasing/Consensus of {2} reads", barcodeName, coarseClusterNumber, readExtents.Length);

            var tpl = new TrialTemplate
            {
                Sequence = consensus,
                StartAdapterBases = 0,
                EndAdapterBases = 0
            };

            // Strip the names from the readExtents for now until I can update the MTC
            // TODO: Update the MultiTemplateConsensus to also use named extents to clean up this code
            var subreadIds = readExtents.Select(n => n.Item1).ToArray();
            var unnamedReadExtents = readExtents
                .Select(n => new Tuple<AlignedSequenceReg, IZmwBases>(n.Item2, n.Item3))
                .ToArray();

            // initial consensus finder - improve the consensus before splitting.
            // The multi-template consensus object tracks the templates we've split, and which reads map to each.
            using (var multiTpl = new MultiTemplateConsensus(unnamedReadExtents, tpl, config))
            {
                multiTpl.ImproveConsensus(3);

                if (DoPhasing && MultiTemplateConsensus.splitThreshold < 10000)
                {
                    // The actual phasing happens here!
                    var nSplits = multiTpl.RecursiveSplit(2);
                    var split = nSplits > 1;

                    if (split)
                    {
                        // FIXME - log number of iterations required and warn if it didn't converge.
                        Log(LogLevel.INFO, "Barcode '{0}', Cluster {1} -- Split into {2} haplotypes", barcodeName, coarseClusterNumber, nSplits);
                    }
                }
                else
                {
                    Log(LogLevel.INFO, "Barcode '{0}', Cluster {1} -- Skipping phasing by user request", barcodeName, coarseClusterNumber, readExtents.Length);
                }

                // Compute the final consensus
                const int maxIterations = 15;
                var iterationsTaken = multiTpl.ImproveConsensus(maxIterations);
                var maxIterationsTaken = iterationsTaken.Aggregate(Math.Max);

                if (maxIterationsTaken < maxIterations)
                {
                    Log(LogLevel.INFO, "Barcode '{0}', Cluster {1} -- Consensus converged in {2} iteration(s)", barcodeName, coarseClusterNumber, maxIterationsTaken);
                }
                else
                {
                    Log(LogLevel.INFO, "Barcode '{0}', Cluster {1} -- Consensus did not converge in {2} iterations", barcodeName, coarseClusterNumber, maxIterationsTaken);
                }

                // Total coverage per template
                var readsPerTpl = multiTpl.MappingPosteriors().Map(v => (int)Math.Round(v.Sum()));
                var stats = multiTpl.ConsensusStats();
                var scorers = sortByCoverage
                    ? multiTpl
                        .Scorers.Select((scorer, idx) => Tuple.Create(scorer, idx))
                        .OrderByDescending(scIdx => readsPerTpl[scIdx.Item2])
                        .Select((scIdx, idx) => Tuple.Create(scIdx.Item1, scIdx.Item2, idx))
                    : multiTpl.Scorers.Select((scorer, idx) => Tuple.Create(scorer, idx, idx));
                var results = scorers.Select(
                scorerIdxPhase =>
                {
                    var scorer = scorerIdxPhase.Item1;
                    var idx = scorerIdxPhase.Item2;
                    var phaseIdx = scorerIdxPhase.Item3;

                    var totalCoverage = readsPerTpl[idx];
                    var tplStats = stats[idx];
                    var mappingScores = multiTpl.MappingPosteriors()[idx];
                    var fastaName = String.Format("Barcode{0}_Cluster{1}_Phase{2}_NumReads{3}", barcodeName, coarseClusterNumber, idx, totalCoverage);

                    var r = new ConsensusResult
                    {
                        CoarseClusterIndex = coarseClusterNumber,
                        PhaseIndex = phaseIdx,
                        BarcodeName = barcodeName,
                        CoarseClusterName = "Cluster" + coarseClusterNumber,
                        FineClusterName = "Phase" + phaseIdx,
                        FastaName = fastaName,
                        SubreadIds = subreadIds,
                        MappingScores = mappingScores, 
                        Sequence = scorer.Template.GetSequence(Strand.Forward),
                        Coverage = totalCoverage,
                        Stats = tplStats,
                        PredictedAccuracy = 1.0f - (float)tplStats.Select(stat => QVs.QVToProb(stat.ConsensusQV)).Average(),
                        ConsensusConverged = iterationsTaken[idx] < maxIterations,
                        ChimeraResult = ChimeraDetector.defaultChimeraHit,
                        DuplicateOf = String.Empty
                    };

                    return r;
                }).ToArray();

                return LabelDuplicates(results);
            }
        }
    }

    /// <summary>
    /// Per-base consensus statistics
    /// </summary>
    public class ConsensusBaseStats
    {
        /// <summary>
        /// Consensus QV of base
        /// </summary>
        public byte ConsensusQV;

        /// <summary>
        /// Base identity
        /// </summary>
        public char Base;

        /// <summary>
        /// Base coverage
        /// </summary>
        public int Coverage;
    }
    
    /// <summary>
    /// Summary of a chimera call
    /// </summary>
    public class ChimeraResult
    {
        public string SequenceA;
        public string SequenceB;
        public int CrossoverPoint;
        public float ChimeraScore;
    }

    /// <summary>
    /// Consensus result data for one haplotype of an amplicon
    /// </summary>
    public class ConsensusResult
    {
        public int CoarseClusterIndex;
        public int PhaseIndex;

        public string BarcodeName;
        public string CoarseClusterName;
        public string FineClusterName;

        public int Coverage;
        public string[] SubreadIds;
        public float[] MappingScores;
        public float PredictedAccuracy;

        public string FastaName;
        public string Sequence;
        public ConsensusBaseStats[] Stats;
        public ChimeraDetector.chimeraHit ChimeraResult;
        public bool ConsensusConverged;
        public string DuplicateOf;
    }
}
