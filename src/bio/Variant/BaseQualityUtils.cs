using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Bio.Variant
{
	public static class BaseQualityUtils
    {   
		static double[] PhredProbCorrectCache;
		static double[] PhredProbIncorrectCache;

		static BaseQualityUtils() {
			PhredProbCorrectCache = Enumerable.Range (0, QualitativeSequence.Phred_MaxQualityScore).
				Select (x => 1 - Math.Pow (10.0, -x / 10.0)).ToArray();
			PhredProbIncorrectCache = PhredProbCorrectCache.Select (x => (1 - x)).ToArray ();
		}

        public static double GetLog10ErrorProbability(int phredScore)
        {
            return phredScore / -10.0;
        }

		/// <summary>
		/// Gets the error probability.
		/// NOT LOGGED.
		/// </summary>
		/// <returns>The error probability.</returns>
		/// <param name="phredScore">Phred score.</param>
		public static double GetErrorProbability(int phredScore)
		{
			return PhredProbIncorrectCache [phredScore];
		}

		/// <summary>
		/// Gets the probability base is correct.
		/// NOT LOGGED.
		/// </summary>
		/// <returns>The error probability.</returns>
		/// <param name="phredScore">Phred score.</param>
		public static double GetCorrectProbability(int phredScore)
		{
			return PhredProbCorrectCache [phredScore];
		}



    }
}
