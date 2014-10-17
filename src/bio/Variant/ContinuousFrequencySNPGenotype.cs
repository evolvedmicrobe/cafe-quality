using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Bio.Variant
{
	public class ContinuousFrequencySNPGenotype : ContinuousFrequencyGenotype
	{
		public BasePairFrequencies Frequencies;

		public int[] OriginalBasePairCounts;

        public int TotalObservedBases { 
            get { return ResultType == GenotypeCallResult.GenotypeCalled ? OriginalBasePairCounts.Sum(): 0; } 
        }

		public string ReferenceName { get; protected set; }

		protected int pOriginalPosition;

		protected int InsertionOffset;

		/// <summary>
		/// Gets the original position of the SNP or Null if an insertion relative to reference.
		/// </summary>
		/// <value>The original position.</value>
		public int? OriginalPosition {
			get{
				if (InsertionOffset == 0)
                    return pOriginalPosition;
				else
                    return null;
					
			}
		}

		#region Constructors
		/// <summary>
		/// Initializes a new instance of the <see cref="Bio.Variant.ContinuousFrequencySNPGenotype"/> class.
		/// </summary>
		/// <param name="res">Res.</param>
		/// <param name="pileup">Pileup.</param>
		public ContinuousFrequencySNPGenotype(GenotypeCallResult res, PileUp pileup = null) {
			ResultType = res;
			if (pileup != null) {
				pOriginalPosition = pileup.Position;
				InsertionOffset = pileup.InsertionOffSet;
			}
		}

		public ContinuousFrequencySNPGenotype(BasePairFrequencies frequencies, int[] originalReadCounts, PileUp pileup = null):
		this(GenotypeCallResult.GenotypeCalled,pileup) {
			Frequencies = frequencies;
			OriginalBasePairCounts = originalReadCounts;
			var freq = frequencies.Frequencies;
			this.indexMax = 0;
			var max = freq [0];
			for (int i = 1; i < freq.Length; i++) {
				if (freq [i] > max) {
					max = freq [i];
					indexMax = i;
				}
			}
		}
		#endregion

		static List<string> basesInOrder = Enumerable.Range(0,4).
			Select(x => BaseAndQuality.Get_DNA_MappingForIndex(x).
				ToString()).
			ToList();


		#region implemented abstract members of ContinuousFrequencyGenotype
		public override List<string> GetGenotypesPresent ()
		{
			return basesInOrder;
		}
		public override List<double> GetGenotypeFrequencies ()
		{
			return Frequencies.Frequencies.ToList ();
		}

		/// <summary>
		/// Keeps track of the index for the highest base.
		/// </summary>
		private int indexMax;

		/// <summary>
		/// Gets the frequency for the most frequent base.
		/// </summary>
		/// <returns>The most frequent genotype.</returns>
		public override string GetMostFrequentGenotype(){
			return Bio.Variant.BaseAndQuality.Get_DNA_MappingForIndex (indexMax).ToString();
		}

		/// <summary>
		/// Gets the highest frequency of any base.
		/// </summary>
		/// <returns>The highest frequency.</returns>
		public override double GetHighestFrequency (){
			return Frequencies.Frequencies [indexMax];
		}
		#endregion
	}
}

