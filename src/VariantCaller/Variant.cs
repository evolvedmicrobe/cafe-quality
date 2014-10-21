using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Bio;

namespace VariantCaller
{
    /// <summary>
    /// Top level class for holding variant information.  
    /// This is implemented in two sub classes, SNPVariant and IndelVariant.
    /// </summary>
    public abstract class Variant
    {
        /// <summary>
        /// What reference genome is this connected to?
        /// </summary>
        public Sequence RefSeq { get; protected set; }

        /// <summary>
        /// SNP, indel, etc.
        /// Now redundant with the Variant Type
        /// </summary>
        public VariantType Type { get; protected set; }

        /// <summary>
        /// 0-based start position of variant.  For Indels, this is the left most position 
        /// BEFORE the event. e.g.
        /// 
        /// AAATTTAAA   -> is a deletion of length 3 starting at position 2.
        /// AAA---AAA
        /// 
        /// A--TTAAA -> is an insertion of length 2 starting at position 0.
        /// ACCTTAAA
        /// </summary>
        public int StartPosition { get; protected set; }

        /// <summary>
        /// O-based end index (same as start for SNPs).
        /// A SNP is of length 1, a deletion of length 2 is 2, etc.
        /// </summary>
        public int Length { get; protected set; }

        /// <summary>
        /// The position in the alignment where this variant ends.
        /// For a SNP, this is the position AFTER the SNP.  Same for indels.
        /// </summary>
        public int EndPosition
        {
            get { return StartPosition + Length; }
        }

        /// <summary>
        /// Is the variant at the very start or end of an alignment? 
        /// (that is it was called based on the first or last base seen on 
        /// either sequence in the alignment.)
        /// These can have certain pathologies so we take note and keep an eye on them.
        /// They should almost always be excluded by the alignment algorithm clipping at the ends.
        /// </summary>
        public bool AtEndOfAlignment { get; protected set; }
        
        public Sequence ReferenceBases  {
            get {
                return (Sequence) RefSeq.GetSubSequence(StartPosition, (EndPosition - StartPosition + 1));
            }
        }

        public Variant(int position, Sequence reference, bool atAlignmentEnd = false)
        {
            StartPosition = position;
            RefSeq = reference;
            AtEndOfAlignment = atAlignmentEnd;
        }
    }
}
