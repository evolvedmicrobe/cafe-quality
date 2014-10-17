using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Bio;
using Bio.Algorithms.MUMmer;
using Bio.Algorithms.Alignment;
using Bio.Extensions;

namespace VariantCaller
{
    public class Reference
    {
        NucmerQueryable nucmer;
        public readonly Sequence ReferenceSequence;
        public Reference(Sequence seq)
        {
            ReferenceSequence = new Sequence(DnaAlphabet.Instance, seq.ToArray(), false);
            nucmer= new NucmerQueryable(ReferenceSequence);


        }
        public IEnumerable<IEnumerable<DeltaAlignment>> GetDeltaAlignments(ISequence querySequence)
        {
            return nucmer.GetDeltaAlignments(querySequence);
            
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
                var seq = ReferenceSequence.GetSubSequence(start, (end - start + 1));
                seq.ID = "Ref:" + start.ToString() + "-" + end.ToString();
                return new SectionOfReferenceGenome() { Start = start, End = end, Seq = seq as Sequence };
            }
            else
            {
                var seq1 = ReferenceSequence.GetSubSequence(start, ReferenceSequence.Count - start);
                var seq2 = ReferenceSequence.GetSubSequence(0, end);
                List<byte> seqs = new List<byte>(seq1);
                seqs.AddRange(seq2);
                var seq = new Sequence(DnaAlphabet.Instance,seqs.ToArray());
                seq.ID = "Ref:" + start.ToString() + "-" + end.ToString();
                return new SectionOfReferenceGenome() {Start = start, End = end, Seq = seq};
            }

        }
        public const string DELTA_ALIGNMENT_METADATAKEY = "DeltaAlns";


		public CompactSAMSequence AlignSequence(ISequence seq)
		{
            return null;
            //CompactSAMSequence css; 
            //if (seq is QualitativeSequence) {
            //    var qs = seq as QualitativeSequence;
            //    css = new CompactSAMSequence (seq.Alphabet,
            //        FastQFormatType.GATK_Recalibrated, 
            //        qs.ToArray (), 
            //        qs.GetPhredQualityScores ());
            //}
            //else
            //{
            //    css = new CompactSAMSequence (seq.Alphabet,
            //        FastQFormatType.GATK_Recalibrated, 
            //        seq.ToArray (), 
            //        Enumerable.Repeat ((byte)30, seq.Count).ToArray ());
            //}

		}

      
       
	
    }
}
