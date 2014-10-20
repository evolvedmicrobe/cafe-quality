using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Bio;

namespace VariantCaller
{
    public class Variant
    {
        /// <summary>
        /// What reference genome is this connected to?
        /// </summary>
        public Reference Genome;

        /// <summary>
        /// SNP, indel, etc.
        /// </summary>
        public VariantType Type;

        /// <summary>
        /// 0-based start position of variant, 
        /// </summary>
        public readonly int StartPosition;

        /// <summary>
        /// O-based end index (same as start for SNPs)
        /// </summary>
        public readonly int EndPosition;

        public Sequence ReferenceBases
        {
            get
            {
                return (Sequence) Genome.ReferenceSequence.GetSubSequence(StartPosition, (EndPosition - StartPosition + 1));
            }
        }


    }
}
