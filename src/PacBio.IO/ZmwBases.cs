using System;
using System.Collections.Generic;
using System.Linq;
using PacBio.Utils;
using PacBio.Align;

namespace PacBio.IO
{
    /// <summary>
    /// ZMW metrics derived from ZmwBases
    /// </summary>
    public class ZmwMetricsBases : IZmwMetricsBases
    {
        private readonly IZmwBases zmwBases;

        public ZmwMetricsBases(IZmwBases b)
        {
            zmwBases = b;

            // Assign default values (not informed) for read-quality assessments
            ReadScore = -1;
            Productivity = ProductivityClass.NotDefined;
            Loading = LoadingClass.NotDefined;

            Regions = new List<RegionAnnotator.Region>();
        }

        /// <summary>
        /// ReadScore should be assigned via some external algorithm
        /// </summary>
        public float ReadScore { get; set; }

        /// <summary>
        /// Productivity should be assigned via some external classifier
        /// </summary>
        public ProductivityClass Productivity { get; set; }

        /// <summary>
        /// LoadingClass is assigned via LoadingClassifier
        /// </summary>
        public LoadingClass Loading { get; set; }

        /// <summary>
        /// ControlReadQual is assigned by aligning to control template.
        /// </summary>
        public float ControlReadQual { get; set; }

        /// <summary>
        /// ControlReadLen is assigned by aligning to control template.
        /// </summary>
        public uint ControlReadLen { get; set; }

        public float[] BaseFraction
        {
            get { return ZmwMetricsBasesFunc.BaseFraction(zmwBases); }
        }

        public float[] CmBasQv
        {
            get { throw new NotImplementedException(); }
        }

        public float RmBasQv
        {
            get { throw new NotImplementedException(); }
        }

        #region Rich QV Metrics

        public float[] CmInsQv
        {
            get { throw new NotImplementedException(); }
        }

        public float[] CmDelQv
        {
            get { throw new NotImplementedException(); }
        }

        public float[] CmSubQv
        {
            get { throw new NotImplementedException(); }
        }

        public float RmInsQv
        {
            get { throw new NotImplementedException(); }
        }

        public float RmDelQv
        {
            get { throw new NotImplementedException(); }
        }

        public float RmSubQv
        {
            get { throw new NotImplementedException(); }
        }

        #endregion

        #region Time-Dependent Metrics

        /// <summary>
        /// Default or Config-defined Time-Dependent Metric interval (sec).
        /// </summary>
        private static float _cfgTdmIntervalSec = 10.0f;
        public static float CfgTdmIntervalSec
        {
            get { return _cfgTdmIntervalSec; }
            set { _cfgTdmIntervalSec = value; }
        }
        
        public float TdmIntervalSec
        {
            get { return CfgTdmIntervalSec; }
        }

        /// <summary>
        /// Default or Config-defined Pause IPD threshold (sec).
        /// </summary>
        private static float _cfgPauseIpdThreshSec = 2.5f;
        public static float CfgPauseIpdThreshSec
        {
            get { return _cfgPauseIpdThreshSec; }
            set { _cfgPauseIpdThreshSec = value; }
        }

        public float PauseIpdThreshSec
        {
            get { return CfgPauseIpdThreshSec; }
        }

        public short[] NumPauseVsT
        {
            get { throw new NotImplementedException(); }
        }

        public short[] NumBaseVsT
        {
            get { throw new NotImplementedException(); }
        }

        public float[] BaseRateVsT
        {
            get { throw new NotImplementedException(); }
        }

        #endregion

        public IList<RegionAnnotator.Region> Regions
        {
            get; private set;
        }


        public float BaseRate
        {
            get { throw new NotImplementedException(); }
        }

        public float BaseWidth
        {
            get { throw new NotImplementedException(); }
        }

        public float BaseIpd
        {
            get { throw new NotImplementedException(); }
        }

        public float Pausiness
        {
            get { throw new NotImplementedException(); }
        }

        public float LocalBaseRate
        {
            get { throw new NotImplementedException(); }
        }

        public float DarkBaseRate
        {
            get { throw new NotImplementedException(); }
        }

        public float[] HQRegionSNR
        {
            get { throw new NotImplementedException(); }
        }

        public float[] HQRegionIntraPulseStd
        {
            get { throw new NotImplementedException(); }
        }

        public float HQRegionStartTime
        {
            get { throw new NotImplementedException(); }
        }

        public float HQRegionEndTime
        {
            get { throw new NotImplementedException(); }
        }

        public float[] HQRegionEstPkmid
        {
            get { throw new NotImplementedException(); }
        }

        public float[] HQRegionEstPkstd
        {
            get { throw new NotImplementedException(); }
        }

        public float[] HQRegionPkzvar
        {
            get { throw new NotImplementedException(); }
        }
    }

    /// <summary>
    /// Provides access to Zmw basecall data. See IZmwBases for field documentatio
    /// </summary>
    public class ZmwBases : IZmwBases
    {
        public uint NumBases
        {
            get
            {
                if (Base == null)
                    return 0;

                return (uint)Base.Count;
            }
        }

        public IList<char> Base
        {
            get;
            set;
        }

        public IList<byte> QV
        {
            get;
            set;
        }

        public IList<int> PulseIndex
        {
            get;
            set;
        }

        public string Sequence
        {
            get { return new String(Base.ToArray()); }
        }

        public bool HaveRichQvData
        {
            get { return (DeletionQV != null); }
        }

        public bool HaveKineticData
        {
            get { return (WidthInFrames != null); }
        }

        public IList<byte> InsertionQV { get; set; }
        public IList<byte> MergeQV { get; set; }

        public IList<byte> DeletionQV { get; set; }
        public IList<char> DeletionTag { get; set; }

        public IList<byte> SubstitutionQV { get; set; }
        public IList<char> SubstitutionTag { get; set; }

        public IList<ushort> WidthInFrames { get; set; }
        public IList<ushort> PreBaseFrames { get; set; }

        public ISequencingZmw Zmw { get; private set; }

        public IZmwMetricsBases Metrics { get; private set; }

        /// <summary>
        /// Use this constructor if you don't have your pulses
        /// </summary>
        /// <param name="zmw"></param>
        public ZmwBases(ISequencingZmw zmw)
        {
            Zmw = zmw;
            Metrics = new ZmwMetricsBases(this);
        }

        /// <summary>
        /// Instantiate the 'Null' object this way if pulses are not available
        /// </summary>
        /// <param name="zmw"></param>
        /// <returns></returns>
        public static ZmwBases Null(ISequencingZmw zmw)
        {
            return new ZmwBases(zmw)
                       {
                           Base = new char[0],
                           PulseIndex = new int[0],
                           QV = new byte[0],
                           InsertionQV = new byte[0],
                           MergeQV = new byte[0],
                           DeletionQV = new byte[0],
                           DeletionTag = new char[0],
                           SubstitutionQV = new byte[0],
                           SubstitutionTag = new char[0],
                           WidthInFrames = new ushort[0],
                           PreBaseFrames = new ushort[0]
                       };
        }
    }

    /// <summary>
    /// Provides access to consensus basecalls data. See IZmwConsensusBases for details
    /// </summary>
    public class ZmwConsensusBases : ZmwBases, IZmwConsensusBases
    {
        public ZmwConsensusBases(ISequencingZmw zmw, IEnumerable<DelimitedSeqReg> regions, float predictedArracy, int estimatedLength, int numAcceptedMutations)
            : base(zmw)
        {
            InsertRegions = regions.ToArray();
            //PredictedAccuracy = predictedArracy;
            this.estimatedLength = estimatedLength;

            // Simple HQ Region covering the CCS read
            var hqReg = new RegionAnnotator.Region
                {
                    Start = 0,
                    End = estimatedLength,
                    HoleNumber = zmw.HoleNumber,
                    Score = (int) (predictedArracy*1000),
                    Type = RegionAnnotator.HqRegionType
                };

            var metrics = Metrics as ZmwMetricsBases;

            metrics.ReadScore = predictedArracy;
            
            metrics.Regions.Add(hqReg);

            NumberOfMutations = numAcceptedMutations;
            
            if(estimatedLength > 5)
                metrics.Productivity = ProductivityClass.Productive;
            else
                metrics.Productivity = ProductivityClass.Empty;
        }

        public DelimitedSeqReg[] InsertRegions
        {
            get;
            private set;
        }

        public int NumberOfMutations { get; set; }

        public int NumPasses
        {
            get { return InsertRegions.Length; }
        }

        public float PredictedAccuracy
        {
            get { return Metrics.ReadScore; }
        }

        private int? estimatedLength;

        public int InsertReadLength
        {
            get
            {
                if (estimatedLength.HasValue)
                    return estimatedLength.Value;

                return (int)NumBases;
            }
        }

        /// <summary>
        /// Instantiate the 'Null' object this way if pulses are not available
        /// </summary>
        /// <param name="zmw"></param>
        /// <returns></returns>
        public static new ZmwConsensusBases Null(ISequencingZmw zmw)
        {
            return new ZmwConsensusBases(zmw, Enumerable.Empty<DelimitedSeqReg>(), 0.0f, 0, -999)
            {
                Base = new char[0],
                PulseIndex = null,
                QV = new byte[0],
                InsertionQV = new byte[0],
                DeletionQV = new byte[0],
                DeletionTag = new char[0],
                SubstitutionQV = new byte[0],
                SubstitutionTag = new char[0],
            };
        }

        public static ZmwConsensusBases MetricsNoSequence(ISequencingZmw zmw, float predictedAccuracy, int estimatedLength)
        {
            return new ZmwConsensusBases(zmw, Enumerable.Empty<DelimitedSeqReg>(), predictedAccuracy, estimatedLength, -999)
            {
                Base = new char[0],
                PulseIndex = null,
                QV = new byte[0],
                InsertionQV = new byte[0],
                DeletionQV = new byte[0],
                DeletionTag = new char[0],
                SubstitutionQV = new byte[0],
                SubstitutionTag = new char[0],
            };
        }

        public static ZmwConsensusBases PassRegion(ISequencingZmw zmw, IZmwBases bases, DelimitedSeqReg reg, float predictedAccuracy, int estimatedLength)
        {
            var r = reg;

            return new ZmwConsensusBases(zmw, new[] { reg }, predictedAccuracy, estimatedLength, -999)
            {
                Base = bases.Base.Slice(r.Start, r.Length),
                PulseIndex = null,
                QV = bases.QV.Slice(r.Start, r.Length),
                InsertionQV = bases.InsertionQV.Slice(r.Start, r.Length),
                DeletionQV = bases.DeletionQV.Slice(r.Start, r.Length),
                DeletionTag = bases.DeletionTag.Slice(r.Start, r.Length),
                SubstitutionQV = bases.SubstitutionQV.Slice(r.Start, r.Length),
                SubstitutionTag = bases.SubstitutionTag.Slice(r.Start, r.Length)
            };
        }
    }
}