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

       
    
     
    }
}