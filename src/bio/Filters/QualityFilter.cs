using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Bio.IO;
using Bio.IO.FastQ;

namespace Bio.Filters
{
    /// <summary>
    /// A class that takes a list of qualitative sequences and trims them until the 
    /// </summary>
    public class FastQZippedQualityTrimmedParser : FastQZippedParser
    {
        private int _trimEndQuality;
        private int _meanRequiredQuality;

        /// <summary>
        /// A class that trims ends and removes low quality reads (usually before assembly)
        /// </summary>
        /// <param name="TrimEndTillOver">Remove all sequences from the end until they are over this QC score</param>
        /// <param name="TrimIfMeanLessThan">Skip the read if the mean quality is less than this value </param>
        public FastQZippedQualityTrimmedParser(string filename, int TrimEndTillOver = 20, int TrimIfMeanLessThan = 20)
            : base(filename)
        {
            _trimEndQuality = 20;
            _meanRequiredQuality = 20;
        }
        public  override IEnumerable<Bio.QualitativeSequence> Parse()
        {
            return FilterSequences(base.Parse());
        }
        /// <summary>
        /// Takes a list of sequences and removes any whose quality is below a certain threshold.
        /// </summary>
        /// <param name="sequencesToFilter">A list of sequences to filter</param>
        /// <returns>Trimmed and Filtered sequences</returns>
        public IEnumerable<QualitativeSequence> FilterSequences(IEnumerable<QualitativeSequence> sequencesToFilter)
        {
            foreach (var seq in sequencesToFilter)
            {
                int[] vals = seq.GetQualityScores();
                int lastAcceptableBase = (int)seq.Count - 1;
                while (lastAcceptableBase > 0)
                {
                    if (vals[lastAcceptableBase] >= _trimEndQuality)
                    {
                        break;
                    }
                    lastAcceptableBase--;
                }
                if (lastAcceptableBase > 0)
                {
                    //check mean
                    double mean = vals.Take(lastAcceptableBase + 1).Average();
                    if (mean > _meanRequiredQuality)
                    {
                        yield return seq.GetSubSequence(0, lastAcceptableBase + 1) as QualitativeSequence;
                    }
                }
            }
        }
    }
}
