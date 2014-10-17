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
        public readonly int Position_Start;

        /// <summary>
        /// O-based end index (same as start for SNPs)
        /// </summary>
        public readonly int Position_End;

        public Sequence Reference
        {
            get
            {
                return (Sequence) Genome.ReferenceSequence.GetSubSequence(Position_Start, (Position_End - Position_Start + 1));
            }
        }


    }
}
