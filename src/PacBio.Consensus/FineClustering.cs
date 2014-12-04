using System;
using System.Collections.Generic;
using System.Linq;
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
