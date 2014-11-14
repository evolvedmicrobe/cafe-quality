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


namespace PacBio.Consensus
{
    using MutationType = ConsensusCore.MutationType;

    /// <summary>
    /// Static methods for generating streams of candidate mutations
    /// </summary>
    public static class GenerateMutations
    {
        /// <summary>
        /// Filter a stream of proposed mutations, only accepting those that fall within some range of another list of mutations
        /// </summary>
        /// <param name="muts">Input mutations to filter</param>
        /// <param name="prev">Previous list of mutations</param>
        /// <param name="spacing">Maximum spacing between an accepted mutation and a mutation in the prev list</param>
        /// <returns>The filtered mutation stream</returns>
        public static IEnumerable<Mutation> PrevMutationFilter(IEnumerable<Mutation> muts, List<Mutation> prev, int spacing)
        {
            // No previous mutations -- we won't find one now either, so bail out
            if(prev.Count == 0)
                yield break;
            
            var prevs = prev.OrderBy(p => p.TemplatePosition).ToArray();

            var shiftedTplPos = new int[prevs.Length];
            var shift = 0;

            // Shift the previous mutations to keep track of their new position in the template
            for (int i = 0; i < prevs.Length; i++)
            {
                var m = prevs[i];
                shiftedTplPos[i] = m.TemplatePosition + shift;

                if (m.Type == MutationType.DELETION)
                    shift--;

                if (m.Type == MutationType.INSERTION)
                    shift++;
            }

            // Only return templates that are in the consider regions
            foreach(var m in muts)
            {
                // Do a binary search to find the closest previous mutation.
                // Let it through if it's within spacing of the current mutation

                var l = 0;
                var r = prevs.Length - 1;

                while (r - l > 1)
                {
                    var mid = (r + l)/2;

                    if (shiftedTplPos[mid] > m.TemplatePosition)
                        r = mid;
                    else
                        l = mid;
                }

                // Find the closest previous mutation
                var dist = Math.Min(Math.Abs(shiftedTplPos[l] - m.TemplatePosition),
                                    Math.Abs(shiftedTplPos[r] - m.TemplatePosition));

                if(dist < spacing)
                    yield return m;
            }
        }
            
        /// <summary>
        /// Enumerate all possible single indel or substitution mutations to the TrialTemplate tpl. The template is 
        /// flanked by presumed correct adapter bases that are not mutated
        /// </summary>
        /// <param name="tpl">The template to generate mutations of</param>
        /// <param name="generateSubstitutions">Generate substitution mutations</param>
        /// <returns>An enumerable of Mutation objects</returns>
        public static IEnumerable<Mutation> GenerateUniqueMutations(TrialTemplate tpl, bool generateSubstitutions = true)
        {
            // Attempt to insert or mismatch every base at every positions.
            // Don't mutate anything in the know template adapter region
            var seq = tpl.GetSequence(Strand.Forward);

            // Start is the first base past the adapter
            var start = Math.Max(1, tpl.StartAdapterBases);

            // End is the fist base of the end adapter
            var end = seq.Length - Math.Max(1, tpl.EndAdapterBases);
            
            for (int i = start; i <= end; i++)
            {
                // homopolyerStart is true if we are on the first base of a hp
                // We will only do a matching insertion or a deletion if we're a 
                // the start of a homopolyer, to prevent retesting the same template
                bool homopolyerStart = !(i > 2 && seq[i-1] == seq[i] && i > start);

                foreach (var b in DNA.Bases)
                {
                    if ((b != seq[i] && b != seq[i-1]) || (homopolyerStart && b != seq[i-1]))
                    {
                        // You are allowed to make an insertion before the first base of the adapter region
                        yield return new Mutation {Base = b, TemplatePosition = i, Type = MutationType.INSERTION};
                    }

                    // Don't generate substitutions if not requested
                    // Don't mutate adapter
                    if (generateSubstitutions && b != seq[i] && i < end)
                        yield return new Mutation { Base = b, TemplatePosition = i, Type = MutationType.SUBSTITUTION };
                }

                // Attempt to delete only the first base of a homopolymer run
                // Don't delete adapter
                if (homopolyerStart && i < end)
                {
                    yield return new Mutation {TemplatePosition = i, Type = MutationType.DELETION, Base = 'A'};
                }
            }
        }
    }
}