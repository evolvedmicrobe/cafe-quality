using System;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Bio.Variant
{
	public class ContinuousFrequencyIndelGenotype : ContinuousFrequencyGenotype
	{
        /// <summary>
        /// Is this an insertion or deletion?
        /// </summary>
        public bool IsInsertion;

        /// <summary>
        /// One based location on rCRS genome (including "N")
        /// </summary>
        public int rCRSPosition;

        public List<string> InsertedBases {get; private set;}

        public List<double> Frequencies { get; private set; }

        public double[] Counts { get; private set; }

        private int MaxIndex;

		public ContinuousFrequencyIndelGenotype (bool insertion, List<string> alleles, double[] counts, int rcrs_position)
		{
            rCRSPosition = rcrs_position;
            IsInsertion = insertion;
            InsertedBases = alleles;
            var total = counts.Sum();
            Frequencies = counts.Select(x => x / total).ToList();
            MaxIndex = 0;
            Counts = counts;
            for (int i = 0; i < Frequencies.Count; i++)
            {
                if (Frequencies[i] > Frequencies[MaxIndex])
                {
                    MaxIndex = i;
                }
            }
		}

        #region implemented abstract members of ContinuousFrequencyGenotype

        public override List<string> GetGenotypesPresent()
		{
            return InsertedBases;
		}

		public override List<double> GetGenotypeFrequencies ()
		{
            return Frequencies;
		}
		public override string GetMostFrequentGenotype(){
            return InsertedBases[MaxIndex];
		}

		public override double GetHighestFrequency (){
            return Frequencies[MaxIndex];
		}

		#endregion
	}
}

