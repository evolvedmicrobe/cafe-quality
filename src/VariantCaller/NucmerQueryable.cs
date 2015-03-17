using System;
using System.Collections.Generic;
using System.Linq;
using Bio;
using Bio.Algorithms.Alignment;
using Bio.Extensions;

namespace VariantCaller
{
	/// <summary>
	/// This class allows one to query a reference genome repeatedly with query sequences, without
	/// having to recreate the suffix array.
	/// 
	/// It is a fork of NucmerPairwiseAligner.cs, with a fair bit of code clean up.  This in itself
	/// is a fork of NucMer (http://mummer.sourceforge.net/manual/).  Nucmer calculates the "delta alignments"
	/// only, and here we calculate that first, then calculate everything else needed to get data (Score, alignment, etc.).
	/// </summary>
    public class NucmerQueryable
    {
        private NUCmer nucmer;
        private Sequence ReferenceSequence;

        public NucmerQueryable(Sequence touse, int lengthOfMUM = 20)
        {
            ReferenceSequence = touse;
            // Note: There are lots of parameters here that could be tuned.
            nucmer= new NUCmer(ReferenceSequence)
           {
                LengthOfMUM = lengthOfMUM,
                MinimumScore = Math.Min((int)(touse.Count * .2), 60)
            };
        }

        public PairwiseAlignedSequence GetLongestAlignment(Sequence toAlign)
        {
            var delts = GetDeltaAlignments (toAlign)
                .SelectMany (x => x).ToList();
            if (delts.Count > 0) {
                var score = new Func<DeltaAlignment, int> ((x) => (int)(x.FirstSequenceEnd - x.FirstSequenceStart + 1));
           
                int max = score (delts [0]);
                var mx = delts [0];
                for (int i = 1; i < delts.Count; i++) {
                    var cur = delts [i];
                    var cur_s = score (cur);
                    if (cur_s > max) {
                        mx = cur;
                        max = cur_s;
                    }
                }
                var aln = NucmerPairwiseAligner.ConvertDeltaToAlignment (mx);
                aln.Score = CalculateScore ((Sequence)aln.FirstSequence, (Sequence)aln.SecondSequence);
                return aln;
            }
            return null;
        }

		public List<PairwiseAlignedSequence> GetAlignments(ISequence toAlign)
		{
            if (Math.Min (toAlign.Count, ReferenceSequence.Count) < nucmer.MinimumScore) {
                var msg = "Bad parameter settings for NucmerPairwiseAligner. " +
                    "Tried to align a reference of length " +ReferenceSequence.Count.ToString() +
                    " to a sequence of length " + toAlign.Count.ToString() +
                    " while requiring a minimum score of MinimumScore = " + nucmer.MinimumScore +
                    ". This will prevent any alignments from being returned.";
                throw new ArgumentException (msg);
            }
			var delts = GetDeltaAlignments (toAlign)
				.SelectMany (x => x);
			var alns = NucmerPairwiseAligner.ConvertDeltasToAlignment (delts).ToList();
			//Now to score them
			alns.ForEach (z => {
                z.Score = CalculateScore ( (Sequence)z.FirstSequence, (Sequence) z.SecondSequence);
			});
			return alns;
		}

        private IEnumerable<IEnumerable<DeltaAlignment>> GetDeltaAlignments(ISequence querySequence)
        {
            IEnumerable<ISequence> querySequences = AddReverseComplementsToSequenceList(querySequence);
            foreach (ISequence qs in querySequences)
            {
                yield return nucmer.GetDeltaAlignments(qs, false, qs.IsMarkedAsReverseComplement());
            }
            
        }

        /// <summary>
        /// Given a list of sequences, create a new list with the original sequence followed
        /// by the Reverse Complement of that sequence.
        /// </summary>
        /// <param name="sequenceList">List of sequence.</param>
        /// <returns>Returns the List of sequence.</returns>
        private IEnumerable<ISequence> AddReverseComplementsToSequenceList(ISequence seq)
        {
            yield return seq;

            ISequence rcSequence = seq.GetReverseComplementedSequence();
            if (rcSequence != null)
            {
                rcSequence.MarkAsReverseComplement();
                yield return rcSequence;
            }

        }
		/// <summary>
		/// Calculate the score of alignment.
		/// </summary>
		/// <param name="referenceSequence">Reference sequence.</param>
		/// <param name="querySequence">Query sequence.</param>
		/// <returns>Score of the alignment.</returns>
		protected int CalculateScore(
			Sequence referenceSequence,
			Sequence querySequence)
		{

			if (referenceSequence == null)
			{
				throw new ArgumentNullException("referenceSequence");
			}

			if (querySequence == null)
			{
				throw new ArgumentNullException("querySequence");
			}

			int index;

			int score = 0;

			// For each pair of symbols (characters) in reference and query sequence
			// 1. If the character are different and not alignment character "-", 
			//      then find the cost from Similarity Matrix
			// 2. If Gap Extension cost needs to be used
			//      a. Find how many gaps exists in appropriate sequence (reference / query)
			//          and calculate the score
			// 3. Add the gap open cost
			for (index = 0; index < referenceSequence.Count; index++)
			{
				byte referenceCharacter = referenceSequence[index];
				byte queryCharacter = querySequence[index];

				if (DnaAlphabet.Instance.Gap != referenceCharacter
					&& DnaAlphabet.Instance.Gap != queryCharacter)
				{
					score += nucmer.SimilarityMatrix[referenceCharacter, queryCharacter];
				}
				else
				{
					int gapCount;
					if (DnaAlphabet.Instance.Gap == referenceCharacter)
					{
						gapCount =  NucmerPairwiseAligner.FindExtensionLength(referenceSequence, index);
					}
					else
					{
						gapCount = NucmerPairwiseAligner.FindExtensionLength(querySequence, index);
					}

					score += nucmer.GapOpenCost + (gapCount * nucmer.GapExtensionCost);

					// move the index pointer to end of extension
					index = index + gapCount - 1;		
				}
			}
			return score;
		}

    }
}
