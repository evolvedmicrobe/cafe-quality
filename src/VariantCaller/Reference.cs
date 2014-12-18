using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Bio;
using Bio.Algorithms.MUMmer;
using Bio.Algorithms.Alignment;
using Bio.Extensions;
using System.Diagnostics;

namespace VariantCaller
{
    [DebuggerDisplay("{RefSeq.ID}")]
    public class Reference
    {
        NucmerQueryable nucmer;
        public readonly Sequence RefSeq;

        public Reference(Sequence seq)
        {
            RefSeq = new Sequence(NoGapDnaAlphabet.Instance, seq.ToArray(), true);
            RefSeq.ID = seq.ID;
            nucmer= new NucmerQueryable(RefSeq, 12);
        }

		public List<PairwiseAlignedSequence> AlignSequence(Sequence toAlign)
		{
            var toR = nucmer.GetAlignments(toAlign);
            toR.Sort ((x, y) => -x.Score.CompareTo (y.Score));
            return toR;
		}

        /// <summary>
        /// Gets a reference sequence based on the 0 indexed values for position in the genom.e
        /// </summary>
        /// <param name="start">0 based index</param>
        /// <param name="end">0 based index, inclusive</param>
        /// <returns></returns>
        public SectionOfReferenceGenome GetReferenceSequenceSection(int start, int end)
        {
            if (end > start)
            {
                var seq = RefSeq.GetSubSequence(start, (end - start + 1));
                seq.ID = "Ref:" + start.ToString() + "-" + end.ToString();
                return new SectionOfReferenceGenome() { Start = start, End = end, Seq = seq as Sequence };
            }
            else
            {
                var seq1 = RefSeq.GetSubSequence(start, RefSeq.Count - start);
                var seq2 = RefSeq.GetSubSequence(0, end);
                List<byte> seqs = new List<byte>(seq1);
                seqs.AddRange(seq2);
                var seq = new Sequence(RefSeq.Alphabet, seqs.ToArray(), false);
                seq.ID = "Ref:" + start.ToString() + "-" + end.ToString();
                return new SectionOfReferenceGenome() {Start = start, End = end, Seq = seq};
            }

        }

        public const string DELTA_ALIGNMENT_METADATAKEY = "DeltaAlns";



      
       
	
    }
}
