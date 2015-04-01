using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Bio;
using Bio.Variant;

namespace VariantCaller
{
    public class SNPVariant : Variant
    {
        /// <summary>
        /// The BP in the reference at this position.
        /// </summary>
        public readonly char RefBP;

        /// <summary>
        /// The variant present at this position.
        /// </summary>
        public readonly char  AltBP;

        /// <summary>
        /// Create a new SNP in the given reference at the given position.
        /// </summary>
        /// <param name="position">0-based position on reference.</param>
        /// <param name="altAllele">The variant present (A, C, G, T)</param>
        /// <param name="reference">The Reference seqeunce.</param>
        public SNPVariant(int position, char altAllele, ISequence reference) :
            base (position, reference)
        {
            AltBP = altAllele;
            Type = VariantType.SNP;
            RefBP = (char) reference[position];
            Length = 1;

        }
    }
}
