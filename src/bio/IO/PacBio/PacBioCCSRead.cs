using System;
using Bio;
using Bio.IO.BAM;
using Bio.IO.SAM;
using System.Collections.Generic;
using System.Linq;

namespace Bio.IO.PacBio
{
    public class PacBioCCSRead
    {
        /// <summary>
        /// A measure of CCS read quality, currently capped at 99.9% (QV30)
        /// </summary>
        public readonly float ReadQuality;

        /// <summary>
        /// The read group.
        /// </summary>
        public readonly string ReadGroup;

        /// <summary>
        /// The number passes.
        /// </summary>
        public readonly int NumPasses;

        /// <summary>
        /// The number of mutations accepted.
        /// </summary>
        public readonly int MutationsAccepted;

        public readonly int MutationsTried;

        private readonly int[] statusCounts;

        /// <summary>
        /// Count of subreads successfully added for consensus generation.
        /// </summary>
        public int ReadsSuccessfullyAdded {
            get{ return statusCounts [0]; }
        }

        /// <summary>
        /// Count of subreads not added to consensus for failing alpha/beta mismatch check.
        /// </summary>
        /// <value>The reads alpha beta mismatch.</value>
        public int ReadsAlphaBetaMismatch {
            get { return statusCounts [1]; }
        }

        /// <summary>
        /// Count of subreads not added to consensus for allocating too much memory.
        /// </summary>
        /// <value>The reads mem fail.</value>
        public int ReadsMemFail {
            get { return statusCounts [2]; }
        }

        /// <summary>
        /// Count of subreads not added to consensus for having too low a Z-score.
        /// </summary>
        /// <value>The reads bad zscore.</value>
        public int ReadsBadZscore {
            get { return statusCounts [3]; }
        }
        public int ReadsOther {
            get{ return statusCounts [4]; }
        }

        public readonly float SnrA;
        public readonly float SnrC;
        public readonly float SnrG;
        public readonly float SnrT;

        /// <summary>
        /// What is the hole number for the ZMW.
        /// </summary>
        public int HoleNumber;

        public float GlobalZscore;

        public float AvgZscore;

        public float ComputingMilliSeconds;

        public readonly float[] ZScores;

        QualitativeSequence Sequence;

        public PacBioCCSRead (SAMAlignedSequence s)
        {
            //TODO: Converting from binary to string and back is beyond silly...

            foreach (var v in s.OptionalFields) {
                if (v.Tag == "sn") {
                    var snrs = v.Value.Split (',').Skip (1).Select (x => Convert.ToSingle (x)).ToArray ();
                    SnrA = snrs [0];
                    SnrC = snrs [1];
                    SnrG = snrs [2];
                    SnrT = snrs [3];
                } else if (v.Tag == "zm") {
                    HoleNumber = (int)Convert.ToInt32 (v.Value);
                } else if (v.Tag == "rq") {
                    ReadQuality = Convert.ToInt32 (v.Value) / 1000.0f;
                } else if (v.Tag == "zg") {
                    GlobalZscore = (float)Convert.ToSingle (v.Value);
                } else if (v.Tag == "za") {
                    AvgZscore = (float)Convert.ToSingle (v.Value);
                } else if (v.Tag == "rs") {
                    statusCounts = v.Value.Split (',').Skip (1).Select (x => Convert.ToInt32 (x)).ToArray ();
                } else if (v.Tag == "np") {
                    NumPasses = Convert.ToInt32 (v.Value);
                } else if (v.Tag == "ms") {
                    ComputingMilliSeconds = Convert.ToSingle (v.Value);                    
                } else if (v.Tag == "mt") {
                    MutationsTried = Convert.ToInt32 (v.Value);   
                } else if (v.Tag == "ma") {
                    MutationsAccepted = Convert.ToInt32 (v.Value);
                } else if (v.Tag == "RG") {
                    ReadGroup = v.Value;
                } else if (v.Tag == "zs") {
                    ZScores = v.Value.Split (',').Skip (1).Select (x => Convert.ToSingle (x)).ToArray ();
                }
            }

            Sequence = s.QuerySequence as QualitativeSequence;

        }
    }
}

