using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using ConsensusCore;
using PacBio.Align;
using PacBio.IO;
using PacBio.IO.Fasta;
using PacBio.Utils;
using PacBio.HDF;
using System.ComponentModel;
using System.Reflection;

namespace PacBio.Consensus
{
    /// <summary>
    /// This class is designed to assign CCS reads to different bins based on SNR and coverage levels.
    /// 
    /// Early analysis has indicated this may greatly improve the analysis.  Right now the decision points
    /// for different bins are entirely based on Nigel eyeballing the data.
    /// </summary>
    [Serializable]
    public class ReadConfigurationAssigner
    {
        /// <summary>
        /// The breakpoints between SNR bins.
        /// </summary>
        public readonly float[] MeanSNRBreakPoints = new float[] {4.0F, 200.0F};

        /// <summary>
        /// Coverage breakpoints
        /// </summary>
        public readonly float[] CoverageBreakPoints = new float[] {3.0F, 20.0F};


        public int NumberOfCoverageGroups { get {return CoverageBreakPoints.Length + 1;}}
        public int NumberOfSNRGroups { get { return MeanSNRBreakPoints.Length + 1; } }

        /// <summary>
        /// The parameters for the assigned group, given by [SNR_GROUP,COVERAGE_GROUP]
        /// </summary>
        QvModelParams[,] parametersForGroups;




        public ReadConfigurationAssigner ()
        {
            parametersForGroups = new QvModelParams[MeanSNRBreakPoints.Length + 1, CoverageBreakPoints.Length + 1];
        }

        public ReadConfigurationAssigner(string fname) : this() {
            
        }

        /// <summary>
        /// Determines if the read is in a specific group, used by training to optimize to determine
        /// if the  some parameters.
        /// </summary>
        /// <returns><c>true</c>, if in group was  read, <c>false</c> otherwise.</returns>
        /// <param name="bases">Bases.</param>
        /// <param name="numberOfInserts">Number of inserts.</param>
        /// <param name="snrGroupIndex">Snr group index.</param>
        /// <param name="coverageGroupIndex">Coverage group index.</param>
        public bool ReadInGroup(IZmwBases bases, int numberOfInserts, int snrGroupIndex, int coverageGroupIndex)
        {
            var meanSNR = bases.Metrics.HQRegionSNR.Average ();
            Console.WriteLine ("Mean SNR: " + meanSNR);
            Console.WriteLine ("Number: " + numberOfInserts);
            var i = assignToSNRGroup (meanSNR);
            var j = assignToCoverageGroup (numberOfInserts);

            return i == snrGroupIndex && j == coverageGroupIndex;
        }

        public QvModelParams GetParameters(IZmwBases bases, int numberOfInserts)
        {
            var meanSNR = bases.Metrics.HQRegionSNR.Average ();
            var i = assignToSNRGroup (meanSNR);
            var j = assignToCoverageGroup (numberOfInserts);
            return parametersForGroups [i, j];
        }

        public PartitionAssignment GetGroupAssignment(float meanSNR, int numberOfInserts)
        {
            var i = assignToSNRGroup (meanSNR);
            var j = assignToCoverageGroup (numberOfInserts);
            return new PartitionAssignment(i, j);
        }

        private int assignToSNRGroup(float meanSNR)
        {
            return assignToGroup (meanSNR, MeanSNRBreakPoints);
        }

        private int assignToCoverageGroup(int numberOfInserts)
        {
            return assignToGroup ((float)numberOfInserts, CoverageBreakPoints);
        }

        private int assignToGroup(float value, float[] breakpoints)
        {
            for (int i = 0; i < breakpoints.Length; i++) {
                if (value < breakpoints [i])
                    return i;
            }
            return breakpoints.Length;
        }

        const string BASE_NAME = "P6-C4";

        public static string GetNameForConfiguration(int snrGroup, int covGroup)
        {
            return BASE_NAME + "_" + snrGroup + "_"+ covGroup;
        }
    }

    public struct PartitionAssignment {
        public readonly int SnrGroup;
        public readonly int CoverageGroup;
        public PartitionAssignment(int snrGroup, int coverageGroup)
        {
            SnrGroup = snrGroup;
            CoverageGroup = coverageGroup;
        }
    }
}

