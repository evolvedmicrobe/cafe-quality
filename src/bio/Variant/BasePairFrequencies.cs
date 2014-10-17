using System;
using System.Linq;

namespace Bio.Variant
{
	public class BasePairFrequencies 
	{
		public const int NUM_BASES = 4;

		/// <summary>
		/// The frequencies of each base, as indexed by the BaseAndQuality structs
		/// indexing scheme for A,C,G,T
		/// </summary>
		public double[] Frequencies = new double[NUM_BASES];

		public double A_Frequency {
			get{ return Frequencies[BaseAndQuality.A_BASE_INDEX]; }
			set { Frequencies [BaseAndQuality.A_BASE_INDEX] = value; }
		}

		public double C_Frequency {
			get{ return Frequencies[BaseAndQuality.C_BASE_INDEX]; }
			set { Frequencies [BaseAndQuality.C_BASE_INDEX] = value; }
		}

		public double G_Frequency {
			get{ return Frequencies[BaseAndQuality.G_BASE_INDEX]; }
			set { Frequencies [BaseAndQuality.G_BASE_INDEX] = value; }
		}

		public double T_Frequency {
			get{ return Frequencies[BaseAndQuality.T_BASE_INDEX]; }
			set { Frequencies [BaseAndQuality.T_BASE_INDEX] = value; }
		}

		public BasePairFrequencies() {
			Frequencies = new double[NUM_BASES];
		}
		/// <summary>
		/// Initializes a new basepair frequencies class with a double array
		/// with values equal to the proportions represented by the BaseAndQuality mapping
		/// of basepairs to values.
		/// </summary>
		/// <param name="frequencies">Frequencies.</param>
		public BasePairFrequencies (double[] frequencies )
		{
			if (frequencies == null || 
				frequencies.Length != NUM_BASES || 
				Math.Abs(frequencies.Sum() - 1.0) > 1e-3) {
				throw new ArgumentException ("frequencies");
			}
			Frequencies = frequencies;
		}
	}
}

