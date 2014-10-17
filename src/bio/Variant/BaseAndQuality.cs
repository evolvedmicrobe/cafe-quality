using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Bio.Variant
{
    /// <summary>
    /// A simple struct that represents a base and quality value from a sequencing read.
    /// </summary>
    public struct BaseAndQuality
    {
        /// <summary>
        /// The base present at this position (either: A, C, G, T , N or '-' encoded as 0-5);
        /// </summary>
        public byte Base;

        /// <summary>
		/// The Phred Score = - 10 * Log_10 (Probability Base is Wrong)
        /// </summary>
		public byte PhredScore;

        /// <summary>
        /// Create a new base/probability pair.
        /// </summary>
        /// <param name="bp">A, C, G, T, N or '-'</param>
		/// <param name="phredScore">The Phred score for the probability that the base is correct. </param>
		public BaseAndQuality(byte bp, byte phredScore)
        {
            var baseIndex = validBases[bp];
			var validData = baseIndex > 0 
				&&  phredScore <=  QualitativeSequence.Phred_MaxQualityScore 
				&&  phredScore >= QualitativeSequence.Phred_MinQualityScore;
            if (!validData)
            {
                throw new ArgumentException("Tried to create a BaseAndQual field with invalid data (BP = " + bp.ToString() + " Qual = " +
					phredScore.ToString() + ").");
            }
            Base = (byte) (baseIndex - 1);
			PhredScore = phredScore;
        }

        /// <summary>
        /// An array with >0 values A, C, G, T, N and '-'
        /// Can be used both to validate if a base is acceptable (element > 0) and
        /// determine a lower order encoding for the base ( = element - 1).
        /// </summary>
        static readonly int[] validBases = new int[byte.MaxValue + 1];

        
		/// <summary>
		/// Initializes the <see cref="Bio.Variant.BaseAndQuality"/> struct.
		/// </summary>
        static BaseAndQuality()
        {
            int i = 1;
			foreach (var bp in BasesInOrder)
            {
                validBases[bp] = i++;
            }
        }

		/// <summary>
		/// A list of valid basepairs as well as their ordering when projected down to a smaller array.
		/// e.g. We change the encoding as follows:
		/// 'A' = 0
		/// 'C' = 1
		/// 'G' = 2
		/// 'T' = 3
		/// 'N' = 4
		/// '-' = 5
		/// Note that these values are stored as constants below as well for direct access.
		/// </summary>
		public static readonly char[] BasesInOrder = new char[] { 'A', 'C', 'G', 'T', 'N', '-' };

		public const byte A_BASE_INDEX = 0;
		public const byte C_BASE_INDEX = 1;
		public const byte G_BASE_INDEX = 2;
		public const byte T_BASE_INDEX = 3;
		public const byte N_BASE_INDEX = 4;
		public const byte GAP_BASE_INDEX = 5;

		/// <summary>
		/// Gets the 0 to 5 index for the nucleotides A, C, G, T, N and '-'
		/// </summary>
		/// <returns>The mapping for nucleotide.</returns>
		/// <param name="bp">Bp.</param>
		public static byte Get_0to5_MappingForNucleotide(byte bp)
		{
			var ind = validBases [bp] - 1; 
			if (ind < 0) {
				throw new ArgumentException ("Byte " + bp.ToString () + " is not A, C, G, T, N or gap character that can be mapped to 0 through 5");
			}
			return (byte)ind;
		}

		public static char Get_DNA_MappingForIndex(int ind)
		{
			switch(ind)
			{
				case A_BASE_INDEX:
					return 'A';
				case C_BASE_INDEX:
					return 'C';
				case G_BASE_INDEX:
					return 'G';
				case T_BASE_INDEX:
					return 'T';
				case GAP_BASE_INDEX:
					return '-';
				case N_BASE_INDEX:
					return 'N';
				default:
					throw new ArgumentException("Index does not map to base!");
			}
		}

    }
}

