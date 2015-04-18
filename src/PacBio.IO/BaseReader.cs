#region Copyright (c) 2010, Pacific Biosciences of California, Inc.
//
// All rights reserved.
//
// THIS SOFTWARE CONSTITUTES AND EMBODIES PACIFIC BIOSCIENCES’ CONFIDENTIAL
// AND PROPRIETARY INFORMATION.
//
// Disclosure, redistribution and use of this software is subject to the
// terms and conditions of the applicable written agreement(s) between you
// and Pacific Biosciences, where “you” refers to you or your company or
// organization, as applicable.  Any other disclosure, redistribution or
// use is prohibited.
//
// THIS SOFTWARE IS PROVIDED BY PACIFIC BIOSCIENCES AND ITS CONTRIBUTORS "AS
// IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
#endregion

using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;
using PacBio.HDF;
using PacBio.IO.ICD;
using PacBio.Utils;
using PacBio.Align;

namespace PacBio.IO
{
    /// <summary>
    /// Implements IZmwMetricsBases, by reading metrics out of BaseReader
    /// </summary>
    internal class HDFZmwMetricsBases : IZmwMetricsBases
    {
        private readonly BaseReader reader;
        private readonly int index;

        public HDFZmwMetricsBases(BaseReader r, int idx)
        {
            reader = r;
            index = idx;
        }

        public float[] BaseFraction
        {
            get { return reader.BaseFraction(index); }
        }

        public float BaseRate
        {
            get { return reader.BaseRate(index); }
        }

        public float BaseWidth
        {
            get { return reader.BaseWidth(index); }
        }

        public float BaseIpd
        {
            get { return reader.BaseIpd(index); }
        }
        public float Pausiness
        {
            get { return reader.Pausiness(index); }
        }

        public float LocalBaseRate
        {
            get { return reader.LocalBaseRate(index); }
        }

        public float DarkBaseRate
        {
            get { return reader.DarkBaseRate(index); }
        }

        public float ReadScore
        {
            get { return reader.ReadScore(index); }
        }

        public float[] HQRegionSNR
        {
            get { return reader.HQRegionSNR(index); }
        }

        public float[] HQRegionIntraPulseStd
        {
            get { return reader.HQRegionIntraPulseStd(index); }
        }

        public float HQRegionStartTime
        {
            get { return reader.HQRegionStartTime(index); }
        }

        public float HQRegionEndTime
        {
            get { return reader.HQRegionEndTime(index); }
        }

        public ProductivityClass Productivity
        {
            get { return reader.Productivity(index); }
        }

        public LoadingClass Loading
        {
            get { return reader.Loading(index); }
        }

        public float ControlReadQual
        {
            get { return reader.ControlReadQual(index); }
        }

        public uint ControlReadLen
        {
            get { return reader.ControlReadLen(index); }
        }

        public float[] CmBasQv
        {
            get { return reader.CmBasQv(index); }
        }

        public float[] CmInsQv
        {
            get { return reader.CmInsQv(index); }
        }

        public float[] CmDelQv
        {
            get { return reader.CmDelQv(index); }
        }

        public float[] CmSubQv
        {
            get { return reader.CmSubQv(index); }
        }

        public float RmBasQv
        {
            get { return reader.RmBasQv(index); }
        }

        public float RmInsQv
        {
            get { return reader.RmInsQv(index); }
        }

        public float RmDelQv
        {
            get { return reader.RmDelQv(index); }
        }

        public float RmSubQv
        {
            get { return reader.RmSubQv(index); }
        }

        public float[] HQRegionEstPkmid
        {
            get { return reader.HQRegionEstPkmid(index); }
        }

        public float[] HQRegionEstPkstd
        {
            get { return reader.HQRegionEstPkstd(index); }
        }

        public float[] HQRegionPkzvar
        {
            get { return reader.HQRegionPkzvar(index); }
        }

        public float TdmIntervalSec
        {
            get { return reader.TdmIntervalSec; }
        }

        public float PauseIpdThreshSec
        {
            get { return reader.PauseIpdThreshSec; }
        }

        public short[] NumPauseVsT
        {
            get { return reader.NumPauseVsT(index); }
        }

        public short[] NumBaseVsT
        {
            get { return reader.NumBaseVsT(index); }
        }

        public float[] BaseRateVsT
        {
            get { return reader.BaseRateVsT(index); }
        }

        public IList<RegionAnnotator.Region> Regions
        {
            get { return reader.Regions(index); }
        }
    }

    /// <summary>
    /// Implements IZmwBases by reading base data out of a base reader.
    /// </summary>
    internal class HDFZmwBases : IZmwBases
    {
        private readonly BaseReader reader;
        private readonly int index;
        private readonly HDFZmwMetricsBases metrics;

        public HDFZmwBases(BaseReader r, int idx)
        {
            reader = r;
            index = idx;
            metrics = new HDFZmwMetricsBases(r, idx);
        }

        public ISequencingZmw Zmw
        {
            get { return reader.ZmwSource[index]; }
        }

        public IZmwMetricsBases Metrics
        {
            get { return metrics; }
        }

        #region Basic Data

        public string Sequence
        {
            get { return new String(Base.ToArray()); }
        }

        public uint NumBases
        {
            get { return reader.ZmwNumEvents[index]; }
        }

        public IList<char> Base
        {
            get { return reader.Base(index); }
        }

        public IList<byte> QV
        {
            get { return reader.QV(index); }
        }

        public IList<int> PulseIndex
        {
            get { return reader.PulseIndex(index); }
        }

        #endregion

        public bool HaveRichQvData
        {
            get { return reader.HaveRichQvData; }
        }

        public IList<byte> InsertionQV
        {
            get { return reader.InsertionQV(index); }
        }

        public IList<byte> DeletionQV
        {
            get { return reader.DeletionQV(index); }
        }

        public IList<char> DeletionTag
        {
            get { return reader.DeletionTag(index); }
        }

        public IList<byte> SubstitutionQV
        {
            get { return reader.SubstitutionQV(index); }
        }

        public IList<byte> MergeQV
        {
            get
            {
                if (reader.HaveMergeQvData)
                    return reader.MergeQV(index);
                
                // Need this satisfied for reducer unit tests.
                // Maybe there is a better solution...
                return new byte[NumBases];
                //return null;
            }
        }

        public IList<char> SubstitutionTag
        {
            get { return reader.SubstitutionTag(index); }
        }

        public bool HaveKineticData
        {
            get { return reader.HaveKineticData; }
        }

        public IList<ushort> WidthInFrames
        {
            get { return reader.WidthInFrames(index); }
        }

        public IList<ushort> PreBaseFrames
        {
            get { return reader.PreBaseFrames(index); }
        }
    }


    internal class HDFZmwConsensusBases : IZmwConsensusBases
    {
        private readonly ConsensusBaseReader reader;
        private readonly int index;

        public HDFZmwConsensusBases(ConsensusBaseReader r, int idx)
        {
            reader = r;
            index = idx;
        }

        public int NumberOfMutations { get { return -999; } }
        public ISequencingZmw Zmw
        {
            get { return reader.ZmwSource[index]; }
        }

        public IZmwMetricsBases Metrics
        {
            get { return null; }
        }

        #region Basic Data

        public string Sequence
        {
            get { return new String(Base.ToArray()); }
        }

        public uint NumBases
        {
            get { return reader.ZmwNumEvents[index]; }
        }

        public IList<char> Base
        {
            get { return reader.Base(index); }
        }

        public IList<byte> QV
        {
            get { return reader.QV(index); }
        }

        public IList<int> PulseIndex
        {
            get { return null; }
        }

        #endregion

        #region Rich QV Extension

        public bool HaveRichQvData
        {
            get { return reader.HaveRichQvData; }
        }

        public IList<byte> InsertionQV
        {
            get { return reader.InsertionQV(index); }
        }

        public IList<byte> MergeQV
        {
            get { return null; }
        }

        public IList<byte> DeletionQV
        {
            get { return reader.DeletionQV(index); }
        }

        public IList<char> DeletionTag
        {
            get { return reader.DeletionTag(index); }
        }

        public IList<byte> SubstitutionQV
        {
            get { return reader.SubstitutionQV(index); }
        }

        public IList<char> SubstitutionTag
        {
            get { return reader.SubstitutionTag(index); }
        }

        #endregion

        public DelimitedSeqReg[] InsertRegions { get; internal set; }

        public int NumPasses
        {
            get { return InsertRegions.Length; }
        }

        
        public float PredictedAccuracy
        {
            get { return 1.0f; } // reader.PredictedAccuracy(index); }
        }

        /*
        public int InsertReadLength
        {
            get { return (int) reader.InsertReadLength(index); }    
        }
        */
        #region IZmwBases Members


        public bool HaveKineticData
        {
            get { return false; }
        }

        public IList<ushort> WidthInFrames
        {
            get { throw new NotImplementedException(); }
        }

        public IList<ushort> PreBaseFrames
        {
            get { throw new NotImplementedException(); }
        }

        #endregion
    }

    /// <summary>
    /// This class provides access to Basecalls stored in PacBio basecalls HDF5 files (*.pls.h5 or *.bas.h5)
    /// Bascalls are exposed through the IBaseSource/IZmwBases interface. See those interfaces for documentation.
    /// </summary>
    [HdfIcdEntry(
            Path = "/PulseData/BaseCalls",
            Detail = "Metrics for all pulses classified as base-incorporation events by the base-caller")
    ]
    public class BaseReader : ReaderBase<IZmwBases>, IBaseSource
    {
        // Functions that convert byte codes read from file to classifier outcomes
        private readonly Func<byte, ProductivityClass> fProductivity;
        private readonly Func<byte, LoadingClass> fLoading;

        #region Interface Control Documentation

        public class Icd : HdfIcd<BaseReader>
        {
            // Define all ICD entries that are not covered by member annotations here
            private static IEnumerable<HdfIcdEntry> Mine()
            {
                return new HdfIcdEntry[]
                           {
                           };
            }

            /// <summary>
            /// Commmon pulse and base ICD entries are specified here, for re-use.
            /// If not a full path, the path will be relative to the class HdfIcdEntry. 
            /// </summary>
            /// <returns></returns>
            public static IEnumerable<HdfIcdEntry> Common()
            {
                return new[]
                {
                    new HdfIcdEntry {Path = "/PulseData", Detail = "Top-level group containing all pulse and base data"},
                    new HdfIcdEntry {Path = "SchemaRevision", Detail = "Version or revision number of the group schema"},
                    new HdfIcdEntry {Path = "Content", Detail = "Content description of the group, as dataset [name, type]"},
                    new HdfIcdEntry {Path = "CountStored", Detail = "The number of ZMW reads contained in each dataset"},
                    new HdfIcdEntry {Path = "QVDecoding", Detail = "Description of the quality-value encoding scheme"},
                    new HdfIcdEntry {Path = "ZMW", Detail = "ZMW bookkeeping markup group"},
                    new HdfIcdEntry {Path = "ZMW/Content", Detail = "Content description of the group, as dataset [name, type]"},
                    new HdfIcdEntry {Path = "ZMWMetrics", Detail = "ZMW metrics markup group"},
                };
            }


            public Icd(bool flatten = false)
                : base(Mine().Concat(Common()), flatten)
            {
                // Units (or LookupTable) for classification metrics have to be assigned dynamically
                HdfIcdEntry val;

                if (TryGetValue("/PulseData/BaseCalls/ZMWMetrics/Productivity", out val))
                    val.Units = MappingForEnum<ProductivityClass>().Aggregate((a,b) => String.Format("{0},{1}",a,b));

                if (TryGetValue("/PulseData/BaseCalls/ZMWMetrics/Loading", out val))
                    val.Units = MappingForEnum<LoadingClass>().Aggregate((a, b) => String.Format("{0},{1}", a, b));
            }
        }

        /// <summary>
        /// Return the ICD in flattened form (for the Writer class).
        /// </summary>
        /// <returns></returns>
        public static ICD.Icd GetIcd() { return new Icd(true); }

        #endregion

        // Data that is acutally stored in the HDF file.
        
        [HdfIcdEntry(Path = "Basecall", Detail = "Called base")]
        public Func<int, char[]> Base;

        [HdfIcdEntry(Detail = "Index into called pulses")]
        public Func<int, int[]> PulseIndex;

        [HdfIcdEntry(Path = "QualityValue", Units = "Phred QV",
                     Detail = "Probability of basecalling error at the current base")]
        public Func<int, byte[]> QV;

        // Optional rich-QV data

        [HdfIcdEntry(Units = "Phred QV", Detail = "Probability of deletion error prior to the current base")]
        public Func<int, byte[]> DeletionQV = ThrowNoRichQv<byte[]>;

        [HdfIcdEntry(Detail = "Likely identity of deleted base (if it exists)")]
        public Func<int, char[]> DeletionTag = ThrowNoRichQv<char[]>;

        [HdfIcdEntry(Units = "Phred QV", Detail = "Probability that the current base is an insertion")]
        public Func<int, byte[]> InsertionQV = ThrowNoRichQv<byte[]>;

        [HdfIcdEntry(Units = "Phred QV", Detail = "Probability of merged-pulse error at the current base")]
        public Func<int, byte[]> MergeQV = ThrowNoRichQv<byte[]>;

        [HdfIcdEntry(Units = "Phred QV", Detail = "Probability of substitution error at the current base")]
        public Func<int, byte[]> SubstitutionQV = ThrowNoRichQv<byte[]>;

        [HdfIcdEntry(Detail = "Most likely alternative base")]
        public Func<int, char[]> SubstitutionTag = ThrowNoRichQv<char[]>;

        // Optional minimal kinetic data

        [HdfIcdEntry(Units = "Frames", Detail = "Duration of the base-incorporation event")]
        public Func<int, ushort[]> WidthInFrames = ThrowNoRichQv<ushort[]>;

        [HdfIcdEntry(Units = "Frames", Detail = "Duration between start of base and end of previous base")]
        public Func<int, ushort[]> PreBaseFrames = ThrowNoRichQv<ushort[]>;

        // Metrics that are storage-optional.  If any of these are requested but do not
        // exist in the file, they will be computed on-the-fly.  Those that depend on 
        // pulse data will fail in the latter case, if IZmwPulses are not provided.
        //
        [HdfIcdEntry(Path = "ZMWMetrics/",
                     Detail = "Base fraction by color channel")]
        public Func<int, float[]> BaseFraction;

        [HdfIcdEntry(Path = "ZMWMetrics/", Units = "1/sec",
                     Detail = "Mean rate of base-incorporation events")]
        public Func<int, float> BaseRate;

        [HdfIcdEntry(Path = "ZMWMetrics/", Units = "sec",
                     Detail = "Mean pulse width of base-incorporation events")]
        public Func<int, float> BaseWidth;

        [HdfIcdEntry(Path = "ZMWMetrics/", Units = "sec",
                     Detail = "Robust estimate of the mean inter-pulse distance (IPD) of base-incorporation events")]
        public Func<int, float> BaseIpd;

        [HdfIcdEntry(Path = "ZMWMetrics/", Detail = "Fraction of pause events over the HQ (sequencing) region")]
        public Func<int, float> Pausiness;

        [HdfIcdEntry(Path = "ZMWMetrics/", Units = "1/sec",
                     Detail = "Robust estimate (excluding pauses) of the mean base-incorporation rate")]
        public Func<int, float> LocalBaseRate;

        [HdfIcdEntry(Path = "ZMWMetrics/", Units = "1/sec",
                     Detail = "Predicted local base rate when the chip is not illuminated")]
        public Func<int, float> DarkBaseRate;

        [HdfIcdEntry(Path = "ZMWMetrics/", Detail = "ZMW productivity classification")]
        public Func<int, ProductivityClass> Productivity;
        private Func<int, byte> productivityRaw;

        [HdfIcdEntry(Path = "ZMWMetrics/", Detail = "ZMW loading classification")]
        public Func<int, LoadingClass> Loading;
        private Func<int, byte> loadingRaw;

        [HdfIcdEntry(Path = "ZMWMetrics/", Detail = "ZMW control read length")]
        public Func<int, uint> ControlReadLen;

        [HdfIcdEntry(Path = "ZMWMetrics/", Detail = "ZMW control read accuracy")]
        public Func<int, float> ControlReadQual;

        [HdfIcdEntry(Path = "ZMWMetrics/", Detail = "Read raw accuracy prediction")]
        public Func<int, float> ReadScore;

        [HdfIcdEntry(Path = "ZMWMetrics/", Detail = "Signal-to-Noise Ratio in the HQ region")]
        public Func<int, float[]> HQRegionSNR;

        [HdfIcdEntry(Path = "ZMWMetrics/", Units = "Counts",
                     Detail = "Standard deviation of intra-pulse signal in the HQ region")]
        public Func<int, float[]> HQRegionIntraPulseStd;

        [HdfIcdEntry(Path = "ZMWMetrics/", Units = "sec",
                     Detail = "Start time of the HQ (sequencing) region")]
        public Func<int, float> HQRegionStartTime;

        [HdfIcdEntry(Path = "ZMWMetrics/", Units = "sec",
                     Detail = "End time of the HQ (sequencing) region")]
        public Func<int, float> HQRegionEndTime;

        [HdfIcdEntry(Path = "ZMWMetrics/", Units = "Counts",
                     Detail = "Average estimated intra-pulse amplitude in the HQ region, prior to pulse-calling")]
        public Func<int, float[]> HQRegionEstPkmid;

        [HdfIcdEntry(Path = "ZMWMetrics/", Units = "Counts",
                     Detail = "Average estimated intra-pulse sigma in the HQ region, prior to pulse-calling")]
        public Func<int, float[]> HQRegionEstPkstd;

        [HdfIcdEntry(Path = "ZMWMetrics/",
                     Detail = "Ratio of (observed intra-pulse variance) / (model-predicted intra-pulse variance) in the HQ region")]
        public Func<int, float[]> HQRegionPkzvar;

        [HdfIcdEntry(Path = "ZMWMetrics/",
                     Detail = "Read-mean (over HQ region if present, otherwise global) of the base QualityValue")] 
        public Func<int, float> RmBasQv;

        [HdfIcdEntry(Path = "ZMWMetrics/", Detail = "Read-mean of the base QualityValue by color channel")]
        public Func<int, float[]> CmBasQv;

        [HdfIcdEntry(Path = "ZMWMetrics/", Detail = "Read-mean of the base InsertionQV by color channel")]
        public Func<int, float[]> CmInsQv = ThrowNoRichQv<float[]>;

        [HdfIcdEntry(Path = "ZMWMetrics/", Detail = "Read-mean of the base DeletionQV by color channel")]
        public Func<int, float[]> CmDelQv = ThrowNoRichQv<float[]>;

        [HdfIcdEntry(Path = "ZMWMetrics/", Detail = "Read-mean of the base SubstitutionQV by color channel")]
        public Func<int, float[]> CmSubQv = ThrowNoRichQv<float[]>;

        [HdfIcdEntry(Path = "ZMWMetrics/", Detail = "Read-mean of the base InsertionQV")]
        public Func<int, float> RmInsQv = ThrowNoRichQv<float>;

        [HdfIcdEntry(Path = "ZMWMetrics/", Detail = "Read-mean of the base DeletionQV")]
        public Func<int, float> RmDelQv = ThrowNoRichQv<float>;

        [HdfIcdEntry(Path = "ZMWMetrics/", Detail = "Read-mean of the base SubstitutionQV")]
        public Func<int, float> RmSubQv = ThrowNoRichQv<float>;

        [HdfIcdEntry(Path = "ZMWMetrics/", Detail = "Number of HQ (sequencing) basecalls by time interval")]
        public Func<int, short[]> NumBaseVsT;

        [HdfIcdEntry(Path = "ZMWMetrics/",
                     Detail = "Number of pause events in HQ (sequencing) region by time interval")]
        public Func<int, short[]> NumPauseVsT;

        [HdfIcdEntry(Path = "ZMWMetrics/", Units = "1/sec",
                     Detail = "Base rate in HQ (sequencing) region by time interval")]
        public Func<int, float[]> BaseRateVsT;

        public Func<int, List<RegionAnnotator.Region>> Regions;

        [HdfIcdEntry(Path = "ZMWMetrics/", Units = "sec",
                     Detail = "Time interval used for Time-Dependent Metrics")]
        public float TdmIntervalSec { get; private set; }

        [HdfIcdEntry(Path = "ZMWMetrics/", Units = "sec",
                     Detail = "Threshold on inter-pulse distance for pause classification")]
        public float PauseIpdThreshSec { get; private set; }

        /// <summary>
        /// Use this property to check whether this reader will serve up rich QV data or no
        /// </summary>
        public bool HaveRichQvData { get; private set; }

        /// <summary>
        /// Use this property to check whether this reader will serve up kinetic data or no
        /// </summary>
        public bool HaveKineticData { get; private set; }

        /// <summary>
        /// Use this property to check whether this reader will serve up an associated
        /// reader or source of consensus bases.
        /// </summary>
        public bool HaveConsensusBases
        {
            get
            {
                var pulseData = (IGroup) Chunks.GetChild("PulseData");

                if (pulseData != null)
                {
                    var consensus = pulseData.GetChild("ConsensusBaseCalls");
                    if (consensus != null)
                        return true;
                }

                return false;
            }
        }

        /// <summary>
        /// Use this property to check whether this reader has direct access to the MergeQV or no
        /// If not, it may still provide the MergeQV via the pulse file, when available, so this
        /// property is not exposed as-is. 
        /// </summary>
        internal bool HaveMergeQvData { get; set; }


        #region Structors

        /// <summary>
        /// Factory method for creating the appropriate type of base source from a bas or pls URI.
        /// Clients should use this method for correct handling of single- or multi-part files.
        /// </summary>
        /// <param name="filename"></param>
        /// <returns></returns>
        public static IBaseSource CreateSource(string filename)
        {
            if (!System.IO.File.Exists(filename))
                throw new FileNotFoundException(String.Format("bas.h5 file not found: {0}", filename), filename);

            var partFileNames = Helpers.GetMultiPartUris(filename);

            if (partFileNames != null)
                return new BaseMultiPartReader(partFileNames);

            return new BaseReader(filename);
        }

        /// <summary>
        /// Create an associated consensus base source if consensus basecalls are available.
        /// Otherwise, return null.
        /// </summary>
        /// <returns></returns>
        public IConsensusBaseSource CreateConsensusBaseSource()
        {
            return HaveConsensusBases
                ? new ConsensusBaseReader(Chunks, FileName)
                : null;
        }

        /// <summary>
        /// Return the sequencing enzyme/chemistry/software-version as a 3-tuple,
        /// usually for the purposes of determining the sequencing chemistry (e.g P4-C2).
        /// Raises ChemistryLookupError in the event of an error or missing metadata.xml
        /// </summary>
        /// <value>The chemistry barcode triple.</value>
        public PacBio.Data.ChemistryBarcodeTriple ChemistryBarcode
        {
            get { return ZmwSource.ChemistryBarcode; }
        }

        /// <summary>
        /// Return the Quiver Sequencing Chemistry string (e.g. P4-C2)
        /// </summary>
        /// <value>The sequencing chemistry.</value>
        public string SequencingChemistry
        {
            get { return ZmwSource.SequencingChemistry; }
        }

        // Set up the methods to return read-level metrics data.
        protected void SetupMetricsLookup()
        {
            // Optional Metrics
            // ReSharper disable RedundantTypeArgumentsOfMethod
            BaseFraction = MakeOptMetricLookupFixedWidth<float>("BaseFraction", null);
            BaseRate = MakeOptMetricLookupSingleton<float>("BaseRate", null);
            BaseWidth = MakeOptMetricLookupSingleton<float>("BaseWidth", null);
            BaseIpd = MakeOptMetricLookupSingleton<float>("BaseIpd", null);
            Pausiness = MakeOptMetricLookupSingleton<float>("Pausiness", null);
            LocalBaseRate = MakeOptMetricLookupSingleton<float>("LocalBaseRate", null);
            DarkBaseRate = MakeOptMetricLookupSingleton<float>("DarkBaseRate", null);

            CmBasQv = MakeOptMetricLookupFixedWidth<float>("CmBasQv", null);
            RmBasQv = MakeOptMetricLookupSingleton<float>("RmBasQv", null);
            HQRegionSNR = MakeOptMetricLookupFixedWidth<float>("HQRegionSNR", null);
            HQRegionIntraPulseStd = MakeOptMetricLookupFixedWidth<float>("HQRegionIntraPulseStd", null);
            HQRegionStartTime = MakeOptMetricLookupSingleton<float>("HQRegionStartTime", null);
            HQRegionEndTime = MakeOptMetricLookupSingleton<float>("HQRegionEndTime", null);
            HQRegionEstPkmid = MakeOptMetricLookupFixedWidth<float>("HQRegionEstPkmid", null);
            HQRegionEstPkstd = MakeOptMetricLookupFixedWidth<float>("HQRegionEstPkstd", null);
            HQRegionPkzvar = MakeOptMetricLookupFixedWidth<float>("HQRegionPkzvar", null);
           
            // These have no fall-back
            ReadScore = MakeOptMetricLookupSingleton("ReadScore", b => -1f);
            
            // An extra level of indirection ensures that we translate byte codes to classifications in a backward-compatible way.
            productivityRaw = MakeOptMetricLookupSingleton<byte>("Productivity", b => (byte) ProductivityClass.NotDefined);
            Productivity = idx => fProductivity(productivityRaw(idx));
            
            loadingRaw = MakeOptMetricLookupSingleton<byte>("Loading", b => (byte) LoadingClass.NotDefined);
            Loading = idx => fLoading(loadingRaw(idx));

            ControlReadLen = MakeOptMetricLookupSingleton<uint>("ControlReadLen", b => 0);
            ControlReadQual = MakeOptMetricLookupSingleton<float>("ControlReadQual", b => 0);
            // ReSharper restore RedundantTypeArgumentsOfMethod

            // Rich Quality Values (Optional)
            if (HaveRichQvData)
            {
                // More optional Metrics
                // ReSharper disable RedundantTypeArgumentsOfMethod
                CmInsQv = MakeOptMetricLookupFixedWidth<float>("CmInsQv", null);
                CmDelQv = MakeOptMetricLookupFixedWidth<float>("CmDelQv", null);
                CmSubQv = MakeOptMetricLookupFixedWidth<float>("CmSubQv", null);

                RmInsQv = MakeOptMetricLookupSingleton<float>("RmInsQv", null);
                RmDelQv = MakeOptMetricLookupSingleton<float>("RmDelQv", null);
                RmSubQv = MakeOptMetricLookupSingleton<float>("RmSubQv", null);
                // ReSharper restore RedundantTypeArgumentsOfMethod
            }

            // Time-Dependent Metrics
            // ReSharper disable RedundantTypeArgumentsOfMethod
            NumBaseVsT = MakeOptMetricLookupFixedWidth<short>("NumBaseVsT", null);
            NumPauseVsT = MakeOptMetricLookupFixedWidth<short>("NumPauseVsT", null);
            BaseRateVsT = MakeOptMetricLookupFixedWidth<float>("BaseRateVsT", null);
            // ReSharper restore RedundantTypeArgumentsOfMethod

            // Lazy Setup for reading the Region table
            regionLists = SetupRegions();
            Regions = i => { return regionLists[i]; };
        }

        // Set up the methods to return base feature data
        protected void SetupFeatureLookup()
        {
            // Standard fields
            var baseRead = MakeFeatureLookup1D<byte>("Basecall");
            Base = index => Encoding.ASCII.GetChars(baseRead(index));

            QV = MakeFeatureLookup1D<byte>("QualityValue");
            PulseIndex = MakeFeatureLookup1D<int>("PulseIndex");
            
            // Rich Quality Values (Optional)
            if (HaveRichQvData)
            {
                InsertionQV = MakeFeatureLookup1D<byte>("InsertionQV");

                DeletionQV = MakeFeatureLookup1D<byte>("DeletionQV");
                var delTagRead = MakeFeatureLookup1D<byte>("DeletionTag");
                DeletionTag = index => Encoding.ASCII.GetChars(delTagRead(index));

                SubstitutionQV = MakeFeatureLookup1D<byte>("SubstitutionQV");
                var subTagRead = MakeFeatureLookup1D<byte>("SubstitutionTag");
                SubstitutionTag = index => Encoding.ASCII.GetChars(subTagRead(index));
            }

            if (HaveKineticData)
            {
                WidthInFrames = MakeFeatureLookup1D<ushort>("WidthInFrames");
                PreBaseFrames = MakeFeatureLookup1D<ushort>("PreBaseFrames");
            }

            // MergeQV (may not exist in some bas.h5 files)
            if (HaveMergeQvData)
            {
                MergeQV = MakeFeatureLookup1D<byte>("MergeQV");
            }
        }

        private List<RegionAnnotator.Region>[] regionLists;

        private List<RegionAnnotator.Region>[] SetupRegions()
        {
            IDataset regionsDataset = null;

            var containerGroup = (IGroup)Chunks.GetChild("PulseData");
             if(containerGroup != null)
                regionsDataset = (IDataset)containerGroup.GetChild("Regions");

            var regLists = ZmwSource.Select((idx, zmw) => new List<RegionAnnotator.Region>()).ToArray();
            
            // If we have a regions dataset, load it, and provide access.
            if (regionsDataset != null)
            {
                var rr = new RegionAnnotator.Reader(containerGroup);
                var reg = rr.Regions;

                foreach(var r in reg)
                {
                    var idx = ZmwSource.GetIndexByHoleNumber(r.HoleNumber);
                    regLists[idx].Add(r);
                }
            }

            return regLists;
        }

        /// <summary>
        /// Construct a BaseReader with a corresponding pulse URI provided explicitly.
        /// Use this constructor when bases are written to a separate (bas.h5) file,
        /// but the corresponding pulse calls (pls.h5) are readily available.
        /// This enables metrics (which may rely on both) to be computed if missing.
        /// External clients should use BaseReader.CreateSource() for support of the multi-part layout as of Release 1.4.
        /// </summary>
        /// <param name="basFile">bas.h5 file</param>
        internal BaseReader(string basFile)
            : base(basFile, "BaseCalls")                
        {

            try
            {
                // Determine capacity for providing rich QV data
                HaveRichQvData = (EventGroup.GetChild("DeletionQV") != null);
                HaveKineticData = (EventGroup.GetChild("WidthInFrames") != null);
                HaveMergeQvData = (EventGroup.GetChild("MergeQV") != null);

                // Provide time-dependent metrics parameters
                TdmIntervalSec = ReadOptMetricAnnotation("TdmIntervalSec", ZmwMetricsBases.CfgTdmIntervalSec);
                PauseIpdThreshSec = ReadOptMetricAnnotation("PauseIpdThreshSec", ZmwMetricsBases.CfgPauseIpdThreshSec);
                
                SetupFeatureLookup();

                // This mapping is true of legacy (pre v1.4) files and should never be changed in this code.
                var legacyProdMapping = new[] {"0:Empty", "1:Productive", "2:Other", "255:NotDefined"};

                // Define the mapping for byte ==> ProductivityClass.
                // If the mapping is not defined in the file, it must be a pre v1.4 Release file.
                // In that case, coerce to the legacy definition.
                //
                fProductivity = DefineClassificationMapping("Productivity", ProductivityClass.NotDefined) ??
                                ByteToEnum(legacyProdMapping, ProductivityClass.NotDefined);

                // Define the mapping for byte ==> LoadingClass.
                // If the mapping is not defined in the file, construct it from the current definition.
                //
                fLoading = DefineClassificationMapping("Loading", LoadingClass.NotDefined) ??
                           ByteToEnum(MappingForEnum<LoadingClass>(), LoadingClass.NotDefined);

                SetupMetricsLookup();

                ZmwEventFunc = idx => () => (IZmwBases)new HDFZmwBases(this, idx);
            }
            catch (IOException e)
            {
                throw new IOException(
                    String.Format("Error Reading HDF Base file: {0}", basFile), e);
            }
        }

        // For use as a delegate getter when no rich QV data is available
        private static TX ThrowNoRichQv<TX>(int i)
        {
            throw new Exception("No Rich QV data available");
        }

        #endregion
    }



    [HdfIcdEntry(
            Path = "/PulseData/ConsensusBaseCalls",
            Detail = "Metrics for circular consensus base calls")
    ]
    public class ConsensusBaseReader : ReaderBase<IZmwConsensusBases>, IConsensusBaseSource
    {
        #region Interface Control Documentation

        public class Icd : HdfIcd<ConsensusBaseReader>
        {
            private static IEnumerable<HdfIcdEntry> Mine()
            {
                return new[]
                {
                    // Manual entries
                    new HdfIcdEntry {Path = "Passes", Detail = "Circular consensus passes data group"},
                    new HdfIcdEntry {Path = "Passes/AdapterHitAfter", Detail = "Flag indicating if an adapter hit was detected at the end of this pass"},
                    new HdfIcdEntry {Path = "Passes/AdapterHitBefore", Detail = "Flag indicating if an adapter hit was detected at the beginning of this pass"},
                    new HdfIcdEntry {Path = "Passes/NumPasses", Detail = "ZMW event-stream counts"},
                    new HdfIcdEntry {Path = "Passes/PassDirection", Detail = "Direction of pass across the SMRTbell"},
                    new HdfIcdEntry {Path = "Passes/PassNumBases", Detail = "Number of bases in circular consensus pass"},
                    new HdfIcdEntry {Path = "Passes/PassStartBase", Detail = "Index of first base in circular consensus pass"},

                    // Handle the BestEstimate reads explicitly, since they use the same ConsensusBaseReader
                    new HdfIcdEntry {Path = "/PulseData/BestEstimateBaseCalls", Detail = "Metrics for the best-estimate base calls of the insert"}
                };
            }

            public Icd(bool flatten = false)
                : base(Mine().Concat(BaseReader.Icd.Common()), flatten)
            {
            }
        }

        /// <summary>
        /// Return the ICD in flattened form (for the Writer class).
        /// </summary>
        /// <returns></returns>
        public static ICD.Icd GetIcd() { return new Icd(true); }

        #endregion

        // Data that is acutally stored in the HDF file.

        [HdfIcdEntry(Path = "Basecall", Detail = "Called base")]
        public Func<int, char[]> Base = ThrowNoFeatureData<char[]>;

        [HdfIcdEntry(Path = "QualityValue", Units = "Phred QV", Detail = "Probability of basecalling error at the current base")]
        public Func<int, byte[]> QV = ThrowNoFeatureData<byte[]>;

        // Optional rich-QV data

        [HdfIcdEntry(Units = "Phred QV", Detail = "Probability of deletion error prior to the current base")]
        public Func<int, byte[]> DeletionQV = ThrowNoRichQv<byte[]>;

        [HdfIcdEntry(Detail = "Likely identity of deleted base (if it exists)")]
        public Func<int, char[]> DeletionTag = ThrowNoRichQv<char[]>;

        [HdfIcdEntry(Units = "Phred QV", Detail = "Probability that the current base is an insertion")]
        public Func<int, byte[]> InsertionQV = ThrowNoRichQv<byte[]>;

        [HdfIcdEntry(Units = "Phred QV", Detail = "Probability of substitution error at the current base")]
        public Func<int, byte[]> SubstitutionQV = ThrowNoRichQv<byte[]>;

        [HdfIcdEntry(Detail = "Most likely alternative base")]
        public Func<int, char[]> SubstitutionTag = ThrowNoRichQv<char[]>;

        [HdfIcdEntry(Path = "ZMWMetrics/", Detail = "Predicted average accuracy of consensus sequence")]
        public Func<int, float> ReadScore;
        
        [HdfIcdEntry(Path = "ZMWMetrics/", Detail = "ZMW productivity")]
        public Func<int, byte> Productivity;

        /// <summary>
        /// Use this property to check whether this reader will serve up (base) feature data
        /// </summary>
        public bool HaveFeatureData { get; private set; }

        /// <summary>
        /// Use this property to check whether this reader will serve up rich QV data or not
        /// </summary>
        public bool HaveRichQvData { get; private set; }

        /// <summary>
        /// Initizalize the reader
        /// </summary>
        protected void SetupFeatureLookup()
        {
            // Standard fields (existence required)
            var baseRead = MakeFeatureLookup1D<byte>("Basecall");
            Base = index => Encoding.ASCII.GetChars(baseRead(index));

            QV = MakeFeatureLookup1D<byte>("QualityValue");

            // Rich Quality Values (Optional)
            if (HaveRichQvData)
            {
                InsertionQV = MakeFeatureLookup1D<byte>("InsertionQV");

                DeletionQV = MakeFeatureLookup1D<byte>("DeletionQV");
                var delTagRead = MakeFeatureLookup1D<byte>("DeletionTag");
                DeletionTag = index => Encoding.ASCII.GetChars(delTagRead(index));

                SubstitutionQV = MakeFeatureLookup1D<byte>("SubstitutionQV");

                var subTagRead = MakeFeatureLookup1D<byte>("SubstitutionTag");
                SubstitutionTag = index => Encoding.ASCII.GetChars(subTagRead(index));
            }
        }

        private DelimitedSeqReg[][] insertRegions;
        private void ReadPassData()
        {
            var passGroup = (IGroup) EventGroup.GetChild("Passes");
            // Pull out the NumPulses data
            var numPassDataset = (IDataset) passGroup.GetChild("NumPasses");
            var numPasses = ((int[]) numPassDataset.Read()).Map(v => v);
            var passStartIndex = numPasses.CumulativeSum();

            var passStartBase = (uint[]) passGroup.ReadDataset("PassStartBase");
            var passNumBases = (uint[]) passGroup.ReadDataset("PassNumBases");
            var passDirection = (byte[]) passGroup.ReadDataset("PassDirection");

            //var adapterHitBefore = (byte[]) passGroup.ReadDataset("AdapterHitBefore");
            //var adapterHitAfter = (byte[]) passGroup.ReadDataset("AdapterHitAfter");

            // Load the insert regions and store them.  They need to get added to the
            // HDFZmwConsensusBases inside the thunk.

            insertRegions = NumZmws.Fill(
                i => {
                    var np = numPasses[i];
                    var jo = passStartIndex[i];

                    var inserts = np.Fill(j => {
                                              var start = (int) passStartBase[jo + j];
                                              var end = (int) (start + passNumBases[jo + j] - 1);
                                              var strand = (Strand) passDirection[jo + j];
                                              return new DelimitedSeqReg(start, end, strand);
                                          });
                    return inserts;
                });
        }

        // For use as a delegate getter when no rich QV data is available
        private static TX ThrowNoRichQv<TX>(int i)
        {
            throw new Exception("No Rich QV data available");
        }

        // For use as a delegate getter when no feature data is available
        private static TX ThrowNoFeatureData<TX>(int i)
        {
            throw new Exception("No base feature data available");    
        }

        #region Structors

        /// <summary>
        /// Factory method for creating the appropriate type of consensus source from a pls/bas URI.
        /// Clients should use this method for correct handling of single- or multi-part files.
        /// </summary>
        /// <param name="filename"></param>
        /// <returns></returns>
        public static IConsensusBaseSource CreateSource(string filename)
        {
            var partUris = Helpers.GetMultiPartUris(filename);

            if (partUris != null)
                return new ConsensusBaseMultiPartReader(partUris);

            return new ConsensusBaseReader(filename);
        }

        internal ConsensusBaseReader(IChunkFile chunks, string filename, string flavor = "ConsensusBaseCalls")
            : base(chunks, filename, flavor)
        {
            Init();
        }

        internal ConsensusBaseReader(string filename, string flavor = "ConsensusBaseCalls")
            : base(filename, flavor)
        {
            Init();
        }

        private void Init()
        {
            try
            {   
                // Determine capacity for prividing (any) base feature data
                HaveFeatureData = (EventGroup.GetChild("Basecall") != null);
                
                // Determine capacity for providing rich QV data
                HaveRichQvData = (EventGroup.GetChild("DeletionQV") != null);

                if (HaveFeatureData)
                {
                    SetupFeatureLookup();
                }

                // Set up the metrics data
                //PredictedAccuracy = MakeOptMetricLookupSingleton("PredictedAccuracy", ZmwMetricsBasesFunc.PredictedAccuracy);
               // InsertReadLength = MakeOptMetricLookupSingleton("InsertReadLength", b => b.NumBases);

                ZmwEventFunc =
                    idx => () =>
                        {
                            var cb = new HDFZmwConsensusBases(this, idx) {InsertRegions = insertRegions[idx]};
                            return (IZmwConsensusBases) cb;
                        };

                // Set up the pass data
                ReadPassData();
            }
            catch (IOException e)
            {
                throw new IOException(
                    String.Format("Error Reading HDF Base file: {0}", FileName), e);
            }
        }

        #endregion       
    }

    /// <summary>
    /// A reader for per-acquisition basecalls stored in PacBio bas.h5 (and/or pls.h5) files, where the 
    /// content is split across multiple files, each containing distinct, contiguous chunks of ZMW hole numbers.  
    /// Basecalls are exposed through the IBaseSource/IZmwBases interface. See those interfaces for documentation.
    /// </summary>
    public class BaseMultiPartReader : MultiPartSource<IZmwBases>, IBaseSource
    {
        public BaseMultiPartReader(IEnumerable<string> uris) :
            base(uris.Select(v => new BaseReader(v) as DataSource<IZmwBases>))
        {
        }

        public BaseMultiPartReader(IEnumerable<BaseReader> readers)
            : base(readers.Select(r => r as DataSource<IZmwBases>).OrderBy(r => r[0].Zmw.HoleNumber))
        {
        }

        public bool HaveConsensusBases
        {
            get
            {
                var part0 = parts[0] as IBaseSource;
                return part0 != null ? part0.HaveConsensusBases : false;
            }
        }

        public IConsensusBaseSource CreateConsensusBaseSource()
        {
            return (HaveConsensusBases
                        ? new ConsensusBaseMultiPartReader(parts.Select(p => p.FileName))
                        : null);
        }

        public PacBio.Data.ChemistryBarcodeTriple ChemistryBarcode
        {
            get
            {
                var part0 = parts[0] as IBaseSource;
                return part0.ChemistryBarcode;
            }
        }

        public string SequencingChemistry
        {
            get
            {
                var part0 = parts[0] as IBaseSource;
                return part0.SequencingChemistry;
            }
        }

        public long GetBaseFileSize()
        {
            var fileSizes = OrderedFileNames.Select(uri =>
                {
                    // No leaky
                    using (var basFile = HDFFile.Open(uri, FileMode.Open, FileAccess.Read))
                    {
                        var baseCalls = (IGroup) basFile.GetChild("/PulseData/BaseCalls");
                        var attribute = baseCalls.GetAttribute("BaseFileSize");
                        return attribute != null ? (long)attribute.Read() : -1;
                    }
                });

            return (fileSizes.All(p => p != -1)) ? fileSizes.Sum() : -1;
        }
    }

    /// <summary>
    /// A reader for per-acquisition basecalls stored in PacBio bas.h5 (and/or pls.h5) files, where the 
    /// content is split across multiple files, each containing distinct, contiguous chunks of ZMW hole numbers.  
    /// Consensus basecalls are exposed through the IConsensusBaseSource / IZmwConsensusBases interface.
    /// See those interfaces for documentation.
    /// </summary>
    public class ConsensusBaseMultiPartReader : MultiPartSource<IZmwConsensusBases>, IConsensusBaseSource
    {
        public ConsensusBaseMultiPartReader(IEnumerable<string> filenames, string flavor = "ConsensusBaseCalls") :
            base(filenames.Select(v => new ConsensusBaseReader(v, flavor) as DataSource<IZmwConsensusBases>))
        {
        }
    }
}
