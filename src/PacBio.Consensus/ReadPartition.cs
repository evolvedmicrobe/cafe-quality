using System;
using System.Collections.Generic;
using System.Linq;
using PacBio.Align;
using PacBio.IO;
using PacBio.Utils;

namespace PacBio.Consensus
{
    /// <summary>
    /// Helper class to keep track of the relevant sequence regions being worked on
    /// </summary>
    public class RegionSets
    {
        public List<DelimitedSeqReg> InsertRegions = new List<DelimitedSeqReg>();
        public List<DelimitedSeqReg> PaddedInsertRegions = new List<DelimitedSeqReg>();
        public List<DelimitedSeqReg> AdapterRegions = new List<DelimitedSeqReg>();

        public void Flip()
        {
            InsertRegions = InsertRegions.Select(r => r.Flip).ToList();
            PaddedInsertRegions = PaddedInsertRegions.Select(r => r.Flip).ToList();
        }
    }

    /// <summary>
    /// Take a read and the SMRTbell adapter sequences and determine the regions of the read corresponding to 
    /// alternating passes through the SMRTbell insert
    /// </summary>
    public class ReadPartition
    {
        public ReadPartition(string adapter)
        {
            Adapter = adapter.ToCharArray();
        }

        public char[] Adapter { get; private set; }

        // What is the minimum insert size to put in the regions table?
        public static int MinInsertSize = 10;

        // What is the minimum accuracy required for a single adapter hit?
        public double MinSingleAdapterHitAccuracy = 0.67;
        
        /// <summary>
        /// Finds a list of SeqReg (template region) structures that represent succesive passes over the insert sequence.
        /// Each pass must be surrounded by complete adapter hits with at least 65% accuracy.
        /// </summary>
        /// <param name="bases">The read basecalls object</param>
        /// <param name="adapterPadBases">Number of adapter bases to pad into the returned read regions</param>
        /// <returns></returns>
        public RegionSets GetPartition(IZmwBases bases, int adapterPadBases)
        {
            // Find the set of all good adapter hits, put in order
            var hits = FindAdapterHits(bases).
                OrderBy(al => al.ReadStartBase).Where(h => h.Accuracy > MinSingleAdapterHitAccuracy).ToArray();

            var hqRegion = bases.HQRegion();

            if (hqRegion == null || hqRegion.Length < MinInsertSize)
                return MakeRegions(new List<SimpleAlignment>(), adapterPadBases);

            // Include adapter hits that touch the HQRegion
            var hqRegionHits = hits.Where(h => h.ReadStartBase + h.ReadLength >= hqRegion.Start && h.ReadStartBase <= hqRegion.End).ToList();

            // just use the hits that were found.
            return MakeRegions(hqRegionHits, adapterPadBases, hqRegion);
        }

        /// <summary>
        /// Find a set of good hits between the complete read and an adapter sequence. 
        /// </summary>
        /// <param name="bases">The read base sequence</param>
        /// <returns>A list of alignments between the base sequence and the adapter</returns>
        public List<SimpleAlignment> FindAdapterHits(IZmwBases bases)
        {
            return FindAdapterHits(bases.Sequence);
        }

        struct AdapterEnd
        {
            public int Pos;
            public double Score;
        }

        internal List<SimpleAlignment> FindAdapterHits(string sequence)
        {
            RowMatrix<AlignMode> path = null;
            RowMatrix<short> m1 = null;
            var minSpacing = (int)(Adapter.Length / 1.1);

            try
            {
                // The scoring penalties and minScore control the sensitivity and fp rate
                // of adapter hits.  Make sure you know what you're doing & run the FalsePositiveTest
                // if you adjust these.
                var penalty = new short[] { -7, -7, -13, 4, -4 };
                var minScore = 0;

                // Fill out an alignment matrix of the adapter vs. the entire base sequence
                m1 = GlobalAlign.GetAlignMatrix(Adapter, sequence, penalty, out path, true, false);

                var lastCol = Adapter.Length;
                // Pull out the scores along the read at the last base of the adapter
                var finalCol = m1.Rows.Fill(i => m1[i, lastCol]);

                // Find a reasonable minimum score of adapter alignments to consider, and pick out those points

                var allEndPoints =
                    finalCol.Select((v, i) => new AdapterEnd { Pos = i, Score = v }).Where(p => p.Score > minScore).ToList();

                // Find the set of best alignment endpoints with some minimum spacing between them
                var endPoints = SpacedSelector.BestItems(allEndPoints, p => p.Pos, p => p.Score, minSpacing);

                // Trace back each of the chosen endpoints, and turn it into a complete alignment
                var alignCells = endPoints.Map(ep => GlobalAlign.RunTraceback(path, ep.Pos, lastCol));
                var alignments =
                    alignCells.Map(
                        ac =>
                        SimpleAlignment.OfAlignCells(ac, new TemplateSpec { Template = new string(Adapter) }, sequence));


                // Return only the alignments above some accuracy))
                return alignments.Where(al => al.Accuracy > 0.66 || FlankingMatchScore(al, sequence) > 10).ToList();
            }
            finally
            {
                // Recycle the row of the alignment matrix to save memory
                if (path != null)
                    path.Dispose();

                if (m1 != null)
                    m1.Dispose();
            }
        }


        public int FlankingMatchScore(SimpleAlignment al1, string sequence)
        {
            var len = 70;

            var l1 = Math.Min(al1.ReadStartBase, len);
            var s1 = sequence.Substring(al1.ReadStartBase - l1, l1);


            var l2 = Math.Min(sequence.Length - (al1.ReadStartBase + al1.ReadLength) - 1, len);

            if (l1 > 10 && l2 > 10)
            {
                // We have a flanking sequence -- use it
                var s2 = sequence.Substring(al1.ReadStartBase + al1.ReadLength, l2);
                s2 = DNA.ReverseComplement(s2);


                var penalty = new short[] { -5, -5, -9, 5 };
                var aa = GlobalAlign.GetGlobalAlignScore(s1, s2, penalty);

                return aa;
            }
            else
            {
                // No flanking sequence -- this adapter runs off the end of the read
                // Just accept this as an adapter hit
                return 100;
            }
        }


        static private RegionSets MakeRegions(List<SimpleAlignment> adapterAlignments, int adapterBases, RegionAnnotator.Region hqRegionTrim = null)
        {
            var regionList = new RegionSets();
        
            regionList.AdapterRegions = adapterAlignments.Select(DelimitedSeqReg.OfAlignment).ToList();
            regionList.InsertRegions = MeasureAllInsertsTrimmed(hqRegionTrim, adapterAlignments, 0);
            regionList.PaddedInsertRegions = MeasureAllInsertsTrimmed(hqRegionTrim, adapterAlignments, adapterBases);

            return regionList;
        }
        
        static List<DelimitedSeqReg> MeasureAllInsertsTrimmed(RegionAnnotator.Region hqRegion, IEnumerable<SimpleAlignment> adapterAlignments, int adapterPadBases)
        {
            var regionList = new List<DelimitedSeqReg>();
            var strand = Strand.Forward;

            // Padded insert start
            Func<SimpleAlignment, int> padStartBase =
                al =>
                {
                    var firstAdapterBase = (al.Template.Template.Length - 1) - adapterPadBases;
                    return al.Cells.Last(c => c.Template == firstAdapterBase).Read + 1;
                };

            // Padded insert end
            Func<SimpleAlignment, int> padEndBase =
                al =>
                {
                    var lastAdapterBase = adapterPadBases;
                    return al.Cells.First(c => c.Template == lastAdapterBase).Read - 1;
                };


            var prevIsAdapter = false;
            SimpleAlignment lastAdapter = null;
            var regStart = hqRegion.Start;
            
            // Add all the regions that terminate in a adapter hit
            foreach (var currentAdapter in adapterAlignments)
            {
                if (currentAdapter.ReadStartBase - regStart > MinInsertSize)
                {
                    if (prevIsAdapter)
                    {
                        var reg = new DelimitedSeqReg(padStartBase(lastAdapter), padEndBase(currentAdapter), strand)
                            {
                                AdapterHitBefore = true,
                                AdapterHitAfter = true
                            };

                        regionList.Add(reg);
                    }
                    else
                    {
                        var reg = new DelimitedSeqReg(regStart, padEndBase(currentAdapter), strand)
                            {
                                AdapterHitBefore = false,
                                AdapterHitAfter = true
                            };
                        regionList.Add(reg);
                    }
                    
                    strand = strand == Strand.Forward ? Strand.Reverse : Strand.Forward;
                }

                lastAdapter = currentAdapter;
                prevIsAdapter = true;
                regStart = currentAdapter.ReadStartBase + currentAdapter.ReadLength;
            }

            // Add the final region
            {
                if (hqRegion.End - regStart > MinInsertSize)
                {
                    if (prevIsAdapter)
                    {
                        var reg = new DelimitedSeqReg(padStartBase(lastAdapter), hqRegion.End, strand)
                            {
                                AdapterHitBefore = true,
                                AdapterHitAfter = false
                            };
                        regionList.Add(reg);
                    }
                    else
                    {
                        var reg = new DelimitedSeqReg(regStart, hqRegion.End, strand)
                            {
                                AdapterHitBefore = false,
                                AdapterHitAfter = false
                            };
                        regionList.Add(reg);
                    }
                }
            }

            return regionList;
        }
    }
}
