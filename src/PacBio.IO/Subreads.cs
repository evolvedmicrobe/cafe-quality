using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using PacBio.Align;
using PacBio.Utils;

namespace PacBio.IO
{
    /// <summary>
    /// Represents a PacBio subread. A sub-sequence of a raw read that is trimmed to the HQRegion, and split on adapter hits.
    /// </summary>
    public class Subread
    {
        /// <summary>
        /// Read that this subread comes from.
        /// </summary>
        public IZmwBases Bases { get; private set; }

        /// <summary>
        /// Region descriptor of this subread
        /// </summary>
        public DelimitedSeqReg Region { get; private set; }

        internal Subread(IZmwBases bases, DelimitedSeqReg region)
        {
            Bases = bases;
            Region = region;
        }

        /// <summary>
        /// Forward strand sequence of this subread
        /// </summary>
        [DebuggerBrowsable(DebuggerBrowsableState.Never)]
        public string FwdSequence
        {
            get
            {
                var seq = Bases.Sequence.Substring(Region.Start, Region.Length);
                return seq;
            }
        }

        /// <summary>
        /// Read score of this subread. This number is a prediction of the raw accuracy of the read
        /// </summary>
        [DebuggerBrowsable(DebuggerBrowsableState.Never)]
        public float ReadScore
        {
            get { return Bases.HQRegion().Score; }
        }

        /// <summary>
        /// SNR of the weakest channel in the read
        /// </summary>
        [DebuggerBrowsable(DebuggerBrowsableState.Never)]
        public float MinSnr
        {
            get { return Bases.Metrics.HQRegionSNR.Min(); }
        }

        /// <summary>
        /// Subread Id string
        /// </summary>
        public string SubreadId
        {
            get
            {
                var id = String.Format("{0}/{1}/{2}_{3}", Bases.Zmw.Movie.MovieName, Bases.Zmw.HoleNumber, Region.Start,
                              Region.End);
                return id;
            }
        }

        public override bool Equals(object obj)
        {
            var other = obj as Subread;

            if (other != null && other.SubreadId == SubreadId)
                return true;

            return false;
        }

        public override int GetHashCode()
        {
            return SubreadId.GetHashCode();
        }

        
        /// <summary>
        /// 
        /// </summary>
        /// <param name="b"></param>
        /// <returns></returns>
        public static DelimitedSeqReg[] SubreadRegions(IZmwBases b)
        {
            var hqRegion = b.HQRegion();

            if (hqRegion.Length < 50)
                return new DelimitedSeqReg[] {};


            var adapterHits = b.Metrics.Regions.
                Where(r => r.Type.Type == "Adapter").
                Where(r => r.End >= hqRegion.Start && r.Start <= hqRegion.End);

            var regionList = new List<DelimitedSeqReg>();
            var strand = Strand.Forward;

            var prevIsAdapter = false;
            RegionAnnotator.Region lastAdapter = null;
            var regStart = hqRegion.Start;

            // Add all the regions that terminate in a adapter hit
            foreach (var currentAdapter in adapterHits)
            {
                if (prevIsAdapter)
                {
                    var reg = new DelimitedSeqReg(lastAdapter.End, currentAdapter.Start, strand)
                    {
                        AdapterHitBefore = true,
                        AdapterHitAfter = true
                    };

                    regionList.Add(reg);
                }
                else
                {
                    var reg = new DelimitedSeqReg(regStart, currentAdapter.Start, strand)
                    {
                        AdapterHitBefore = false,
                        AdapterHitAfter = true
                    };
                    regionList.Add(reg);
                }

                strand = strand == Strand.Forward ? Strand.Reverse : Strand.Forward;


                lastAdapter = currentAdapter;
                prevIsAdapter = true;
                regStart = currentAdapter.Start + currentAdapter.Length;
            }

            // Add the final region

            if (prevIsAdapter)
            {
                var reg = new DelimitedSeqReg(lastAdapter.End, hqRegion.End, strand)
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

            return regionList.ToArray();
        }
    }


    /// <summary>
    /// A Subread corresponding to an aligned region
    /// </summary>
    public class AlignedSubread
    {
        public IAlnSummary Alignment { get; private set; }

        public IZmwBases Bases { get; private set; }
        public DelimitedSeqReg Region { get; private set; }

        internal AlignedSubread(IZmwBases bases, DelimitedSeqReg region, IAlnSummary alignment)
        {
            Bases = bases;
            Region = region;
            Alignment = alignment;
        }

        public string FwdSequence
        {
            get
            {
                var seq = Bases.Sequence.Substring(Region.Start, Region.Length);

                if (Alignment.Strand == Strand.Reverse)
                    seq = DNA.ReverseComplement(seq);

                return seq;
            }
        }
    }

    /// <summary>
    /// Helper class
    /// </summary>
    public static class BasesHelpers
    {
        /// <summary>
        /// Extension method for computing the subreads from a read
        /// </summary>
        public static Subread[] Subreads(this IZmwBases bases)
        {
            var reg = Subread.SubreadRegions(bases);
            return reg.Map(r => new Subread(bases, r));
        }
    }
}
