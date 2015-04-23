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
using System.Linq;
using System.Text.RegularExpressions;
using PacBio.Align;

namespace PacBio.IO
{
    /// <summary>
    /// A structure for uniquely identifying a single read. Provides various methods for converting to other
    /// string representations of the same data (eg. a FASTA read tag or trace Url).
    /// There is a difference between SF and Astro read Ids. See: http://smrtwiki/wiki/SpringfieldReadIDProposal
    /// for details.
    /// </summary>
    public struct TraceReference
    {
        public int X;
        public int Y;
        public string RunCode;
        public string MovieFile;

        private static readonly Regex r;

        static TraceReference()
        {
            r = new Regex(@"x(\d+)_y(\d+)_(\d+-\d+)_([^.|]*)", RegexOptions.Compiled);
        }
        
        public static TraceReference Parse(string readTag)
        {
            //x13_y57_1220215-0011_m081218_201142_Uni_p1_b25

            var matches = r.Match(readTag);

            var tref = new TraceReference
                           {
                               X = Int32.Parse(matches.Groups[1].Value),
                               Y = Int32.Parse(matches.Groups[2].Value),
                               RunCode = matches.Groups[3].Value,
                               MovieFile = matches.Groups[4].Value
                           };

            return tref;
        }

        // Springfield FASTA read id: See http://smrtwiki/wiki/SpringfieldReadIDProposal
        public static string CreateSpringfield(IZmwBases read)
        {
            return String.Format("{0}/{1}", read.Zmw.Movie.MovieName, read.Zmw.HoleNumber);
        }

        // Springfield FASTA read id: See http://smrtwiki/wiki/SpringfieldReadIDProposal
        public static string CreateSpringfieldCCS(IZmwBases read)
        {
            return String.Format("{0}/{1}/ccs", read.Zmw.Movie.MovieName, read.Zmw.HoleNumber);
        }

        // Springfield FASTA subread id, with RQ annotation.
        public static string CreateSpringfieldSubread(IZmwBases read, int start, int end)
        {
            var rq = read.Metrics.ReadScore;
            var tagBase = CreateSpringfield(read);

            return String.Format("{0}/{1}_{2} RQ={3:F3}", tagBase, start, end, rq);
        }
    }

    /// <summary>
    /// A structure for uniquely identifying a single read. Provides various methods for converting to other
    /// string representations of the same data (eg. a FASTA read tag or trace Url).
    /// There is a difference between SF and Astro read Ids. See: http://smrtwiki/wiki/SpringfieldReadIDProposal
    /// for details.
    /// </summary>
    public struct SpringfieldTraceReference
    {
        public int HoleNumber;
        public string MovieFile;

        public int SubreadStart;
        public int SubreadEnd;
        public float ReadScore;

        private static readonly Regex r;

        static SpringfieldTraceReference()
        {
            //r = new Regex(@"([^/]*)/([0-9]+)", RegexOptions.Compiled);

            // Bug 19397: Make this compatible with subreads FASTA
            r = new Regex(@"([^/]*)/([0-9]+)(/([0-9]+)_([0-9]+))?( RQ=([0-9]*\.[0-9]+))?", RegexOptions.Compiled);
        }

        public static SpringfieldTraceReference Parse(string readTag)
        {
            var m = r.Match(readTag);

            if (m.Success)
            {
                //for (int i = 3; i < m.Groups.Count; ++i)
                //    Console.WriteLine("Group{0}={1}", i, m.Groups[i].Value);

                var tref = new SpringfieldTraceReference
                               {
                                   MovieFile = m.Groups[1].Value,
                                   HoleNumber = Int32.Parse(m.Groups[2].Value),

                                   // Subread range
                                   SubreadStart = m.Groups[3].Success ? Int32.Parse(m.Groups[4].Value) : 0,
                                   SubreadEnd = m.Groups[3].Success ? Int32.Parse(m.Groups[5].Value) : 0,
                                   
                                   // Read score
                                   ReadScore = m.Groups[6].Success ? Single.Parse(m.Groups[7].Value) : 0,
                               };
                return tref;
            }

            throw new ApplicationException("Failed to match Springfield FASTA header");            
        }

        // Springfield FASTA read id: See http://smrtwiki/wiki/SpringfieldReadIDProposal
        public static string Create(ISequencingZmw zmw)
        {
            return String.Format("{0}/{1}", zmw.Movie.MovieName, zmw.HoleNumber);
        }
    }


    /// <summary>
    /// A class for identifying the type of ZMW being worked with.  The trace file tags each Zmw with a type indicating
    /// whether it is a normal sequencing ZMW, one of various types of fiducial ZMWs, or contains bad data for some reason.
    /// Generally you can use the static members of this class such 
    /// </summary>
    public class ZmwType
    {
        /// <summary>
        /// Name of ZmwType
        /// </summary>
        public string Name { get; private set; }

        /// <summary>
        /// Zmw type code
        /// </summary>
        public byte Code { get; private set; }

        public override string ToString()
        {
            return Name;
        }

        public ZmwType(byte code, string name)
        {
            Code = code;
            Name = name;
        }

        public static ZmwType[] Types;

        static ZmwType()
        {
            var types = new[]
                            {
                                "SEQUENCING",   // Standard sequencing ZMW
                                "ANTIHOLE",     // No ZMW present, used for measuring backgrounds
                                "FIDUCIAL",     // Non standard ZMW - possibly larger than normal
                                "SUSPECT",      // No valid data generated by this ZMW
                                "ANTIMIRROR",   // No micromirror
                                "FDZMW",        // Extra large Zmw
                                "FBZMW",        // FB Zmw?
                                "ANTIBEAMLET",  // No Laser hitting ZMW
                                "OUTSIDEFOV"    // Outside field of view
                            };
            SetupTypes(types);
        }

        public static void SetupTypes(string[] zmwTypeNames)
        {
            Types = zmwTypeNames.Select((name, idx) => new ZmwType((byte)idx, name)).ToArray();
        }

        /// <summary>
        /// A standard sequencing ZMW
        /// </summary>
        public static ZmwType Sequencing
        {
            get { return Types[0]; }
        }

        /// <summary>
        /// No ZMW present at this location, used for measuring backgrounds
        /// </summary>
        public static ZmwType AntiHole
        {
            get { return Types[1]; }
        }

        /// <summary>
        /// A non standard ZMW location may not contain a ZMW
        /// </summary>
        public static ZmwType Fiducial
        {
            get { return Types[2]; }
        }

        /// <summary>
        /// No valid data generated ffor this ZMW. Usually due to interference with camera dead zones
        /// </summary>
        public static ZmwType InvalidData
        {
            get { return Types[3]; }
        }

        /// <summary>
        /// A ZMW without a micromirror
        /// </summary>
        public static ZmwType AntiMirror
        {
            get { return Types[4]; }
        }

        /// <summary>
        /// A large ZMW with an elevated diffusion background signal
        /// </summary>
        public static ZmwType FDZmw
        {
            get { return Types[5]; }
        }

        /// <summary>
        /// FBZmw?
        /// </summary>
        public static ZmwType FBZmw
        {
            get { return Types[6]; }
        }

        /// <summary>
        /// No laser illumination on this ZMW
        /// </summary>
        public static ZmwType AntiBeamlet
        {
            get { return Types[7]; }
        }

        /// <summary>
        /// Ouside high-quality field of view
        /// </summary>
        public static ZmwType OutsideFOV
        {
            get { return Types[8]; }
        }

        public static explicit operator ZmwType(byte code)
        {
            if (!(code < Types.Length))
                throw new ArgumentException(
                    String.Format("Specificed ZmwType: {0} is not registered -- check HoleStatus table", code),
                    "code");

            return Types[code];
        }

        public static explicit operator ZmwType(int code)
        {
            if (!(code < Types.Length))
                throw new ArgumentException(
                    String.Format("Specificed ZmwType: {0} is not registered -- check HoleStatus table", code),
                    "code");

            return Types[code];
        }

        public static explicit operator byte(ZmwType zmwType)
        {
            return zmwType.Code;
        }

        // Comparison operators for ZmwType objects
        public static bool operator ==(ZmwType v1, ZmwType v2)
        {
            return v1.Code == v2.Code;
        }

        public static bool operator !=(ZmwType v1, ZmwType v2)
        {
            return v1.Code != v2.Code;
        }

        public static bool operator ==(ZmwType v1, int v2)
        {
            return v1.Code == v2;
        }

        public static bool operator !=(ZmwType v1, int v2)
        {
            return v1.Code != v2;
        }

        public override bool Equals(object obj)
        {
            var zt = obj as ZmwType;

            if(zt != null && zt.Code == Code)
                return true;

            return false;
        }

        public override int GetHashCode()
        {
            return Code;
        }
    }

    /// <summary>
    /// This will capture the basic context of a Sequencing Trace.  I think it should be defined as all the data
    /// you'd need to make sense of trace and pulse data (ie FrameRate, Duration and Basemap) and a few more basic
    /// items like background and ZMW location. ZMW locations in pixels?
    /// </summary>
    public interface ISequencingZmw
    {
        /// <summary>
        /// The Hole number of the ZMW - this number is a 0 based index into a 2d row major grid of ZMWs.
        /// On Astro this grid was 93x33.
        /// </summary>
        int HoleNumber { get; }

        /// <summary>
        /// The 1-based X position (column number) of the ZMW
        /// </summary>
        int X { get; }

        /// <summary>
        /// The 1-based Y position (row number) of the ZMW
        /// </summary>
        int Y { get; }

        /// <summary>
        /// The look number of the ZMW
        /// </summary>
        Int16 HoleChipLook { get; }

        /// <summary>
        /// Indicates what type of data is available from this ZMW
        /// </summary>
        ZmwType ZmwType { get; }
        /// <summary>
        /// A pointer to the metadata about the movie containing this ZMW
        /// </summary>
        IMovieMetadata Movie { get; }
    }

    /// <summary>
    /// A simple descriptor of an fluorescent analog used in a sequencing experiment.
    /// This type can used to easily derive wavelength maps, and basemaps etc.
    /// </summary>
    public interface IAnalogSpec
    {
        /// <summary>
        /// Identity of the base
        /// </summary>
        char Base { get; }
        /// <summary>
        /// Canonical dye wavelength moniker
        /// </summary>
        float Wavelength { get; }
        /// <summary>
        /// String describing fluorescent moiety
        /// </summary>
        string Label { get; }
        /// <summary>
        /// String describing analog molecule
        /// </summary>
        string Nucleotide { get; }
        /// <summary>
        /// Analog concentration in Moles/L
        /// </summary>
        double Concentration { get; }
    }

    /// <summary>
    /// All the relevant information about a series of pulse calls for a trace.
    /// IList is used because it provides an Item indexer, and we can implement this interface
    /// with wrapper objects that compute derived quantities (i.e. PrePulseTime, or frame -> time conversions).
    /// Need to decide on units, and any associated naming conventions
    /// For signal quantities I would advocate detected photons as the base unit.
    /// To save on space I'm suggesting we use 4-byte floating point numbers instead of doubles.
    /// This information will only be available after a pulse calling stage
    /// </summary>
    public interface IZmwPulses : IBasicZmwPulses
    {
        /// <summary>
        /// Base identity of the observed pulse.  ASCII code representing the observed base.
        /// </summary>
        IList<char> Base { get; }
        
        /// <summary>
        ///  Time in movie when pulse began
        /// </summary>
        IList<float> StartTime { get; }

        /// <summary>
        /// Duration of pulse in seconds
        /// </summary>
        IList<float> WidthInSeconds { get; }

        /// <summary>
        /// Amount of time between end of previous pulse and beginning of current pulse
        /// </summary>
        IList<float> PrePulseTime { get; }

        /// <summary>
        /// Amount of time between end of current pulse and beginning of next pulse
        /// </summary>
        IList<float> PostPulseTime { get; }

        /// <summary>
        /// Overall Phred-scale QV from the pulse classifier
        /// </summary>
        IList<byte> ClassifierQV { get; }

        IList<double> SNR { get; }
        IList<double> PNR { get; }

        IList<float> MinChi2 { get; }
        IList<float> DeltaMinChi2 { get; }
        IList<float> LocalPulseDensity(double timeWindow);

        Func<int, int, IList<float>> PkmidEstProvider { get; }
        Func<int, int, IList<float>> PkstdEstProvider { get; }
        Func<int, float, float> ModelPredictedIntraPulseVariance { get; }
    }

    /// <summary>
    /// All the basic information about a series of pulse calls for a trace.
    /// This interface includes only the data that is directly stored in the pulse file.
    /// The IZmwPulses derived interface includes a bunch of derived metrics.
    /// IList is used because it provides an Item indexer, and we can implement this interface
    /// with wrapper objects that compute derived quantities (i.e. PrePulseTime, or frame -> time conversions).
    /// Need to decide on units, and any associated naming conventions
    /// For signal quantities I would advocate detected photons as the base unit.
    /// To save on space I'm suggesting we use 4-byte floating point numbers instead of doubles.
    /// This information will only be available after a pulse calling stage
    /// </summary>
    public interface IBasicZmwPulses : IZmwEvents
    {
        /// <summary>
        /// Number of pulses in the trace
        /// </summary>
        uint NumPulses { get; }

        /// <summary>
        /// Channel of the observed pulse. Channels are number starting at 0, ordered by ascending emission wavelength
        /// </summary>
        IList<byte> Channel { get; }

        /// <summary>
        /// Frame of movie when pulse began
        /// </summary>
        IList<uint> StartFrame { get; }

        /// <summary>
        /// Duration of pulse in frames
        /// </summary>
        IList<UInt16> WidthInFrames { get; }

        /// <summary>
        /// Phred scale quality value representing the probability of a miscall,
        /// including an insertion error.  Pulses included from an initial no-pulse
        /// label are assigned a value that represents the estimated probability
        /// of deletion error -- for the NO PULSE labeling.
        /// </summary>
        IList<byte> LabelQV { get; }

        /// <summary>
        /// Phred scale quality value representing the probability that a merge
        /// error was made on this pulse, based on the original classification.
        /// Pulses included from an initial no-pulse label are assigned a value of 0.
        /// </summary>
        IList<byte> MergeQV { get; }

        /// <summary>
        /// Did the pulse classifier label the region as a pulse
        /// </summary>
        IList<bool> IsPulse { get; }

        /// <summary>
        /// The maximum signal over the frames of the pulse in each channel. Has dimensions
        /// of NumPulses by NumChannels. The units are always in photons.
        /// </summary>
        IList<UInt16> MaxSignal { get; }

        /// <summary>
        /// The mean signal over the frames of the pulse in each channel. Has dimensions
        /// of NumPulses by NumChannels. The units are always in photons.
        /// </summary>
        UInt16[,] MeanSignal { get; }

        /// <summary>
        /// The mean signal over the interior frames of the pulse in each channel. The interior
        /// frames means that the first and last frames of the the pulse are excluded. If the pulse
        /// has 1 or 2 frames then MidSignal is zero. Has dimensions of NumPulses by NumChannels.
        /// The units are always in photons.
        /// </summary>
        IList<UInt16> MidSignal { get; }

        /// <summary>
        /// The Chi-Squared of the total light of each pulse against the spectral template of each template.
        /// Has dimensions of NumPulses by NumChannels.
        /// </summary>
        float[,] Chi2 { get; }

        /// <summary>
        /// The standard deviation of the signal level in the interior frames of the pulse
        /// </summary>
        IList<UInt16> MidStdDev { get; }

        /// <summary>
        /// The mean bias of the baseline signal, as estimated by the trace signal processing algorithm.
        /// The units are always in photons.
        /// </summary>
        double[] BaselineBias { get; }

        /// <summary>
        /// The standard deviation of the baseline signal (i.e., the average noise level),
        /// as estimated by the trace signal processing algorithm.
        /// The units are always in photons.
        /// </summary>
        double[] BaselineSigma { get; }

        /// <summary>
        /// ClippedFraction (low-end) per channel
        /// </summary>
        float[] ClippedFractionLow { get; }

        /// <summary>
        /// ClippedFraction (high-end) per channel
        /// </summary>
        float[] ClippedFractionHigh { get; }

        /// <summary>
        /// Whether ZMW failed clipped fraction low.
        /// </summary>
        bool FailedClippedFractionLow { get; }

        /// <summary>
        /// Whether ZMW failed clipped fraction high.
        /// </summary>
        bool FailedClippedFractionHigh { get; }

        /// <summary>
        /// The estimate of pulse signal amplitude (pkmid), made as input to pulse classification. 
        /// </summary>
        float[] SignalLevel { get; }

        /// <summary>
        /// The estimate of pulse signal standard deviation (pkstd), made as input to pulse classification.
        /// </summary>
        float[] SignalSigma { get; }

        /// <summary>
        /// A slot for development-related checksum output
        /// </summary>
        double[] Checksum { get; }

        /// <summary>
        /// A reference to the IZmwMetricsPulses object that provides pulse-related summary metrics.
        /// </summary>
        IZmwMetricsPulses Metrics { get; }
    }

    /// <summary>
    /// Pulse-related summary metrics.
    /// These metrics are "storage optional" -- providers that are backed by persistent storage
    /// should handle these properties with methods that fall back to on-the-fly computation
    /// when the data are not present.
    /// </summary>
    public interface IZmwMetricsPulses
    {
        /// <summary>
        /// The mean signal-to-noise ratio by channel
        /// </summary>
        float[] Snr { get; }

        /// <summary>
        /// The average (global) pulse rate in the ZMW (pulses/sec)
        /// </summary>
        float PulseRate { get; }

        /// <summary>
        /// The mean pulse width in the ZMW (sec)
        /// </summary>
        float PulseWidth { get; }

        /// <summary>
        /// Fraction of frames clipped low-end by channel
        /// </summary>
        float[] ClippedFractionLow { get; }

        /// <summary>
        /// Fraction of frames clipped high-end by channel
        /// </summary>
        float[] ClippedFractionHigh { get; }

        /// <summary>
        /// Whether ZMW failed clipped fraction low.
        /// </summary>
        bool FailedClippedFractionLow { get; }

        /// <summary>
        /// Whether ZMW failed clipped fraction high.
        /// </summary>
        bool FailedClippedFractionHigh { get; }
    }

    /// <summary>
    /// All the relevant information about simulated pulses that were used to generate a simulated trace
    /// </summary>
    public interface IZmwPulsesSim : IBasicZmwPulsesSim
    {
        /// <summary>
        /// Base identity of the observed pulse.  ASCII code representing the observed base.
        /// </summary>
        IList<char> Base { get; }

        /// <summary>
        ///  Time in movie when pulse began
        /// </summary>
        IList<double> StartTime { get; }

        /// <summary>
        /// Duration of pulse in seconds
        /// </summary>
        IList<double> WidthInSeconds { get; }

        /// <summary>
        /// Amount of time between end of previous pulse and beginning of current pulse
        /// </summary>
        IList<double> PrePulseTime { get; }

        /// <summary>
        /// Amount of time between end of current pulse and beginning of next pulse
        /// </summary>
        IList<double> PostPulseTime { get; }
    }

    /// <summary>
    /// All the ground-truth information about simulated pulses.
    /// This interface includes only the data that is directly stored in the groung-truth
    /// pulse group, which is typically included in a simulated trace file.
    /// </summary>
    public interface IBasicZmwPulsesSim
    {
        /// <summary>
        /// Number of pulses in the trace
        /// </summary>
        uint NumPulses { get; }

        /// <summary>
        /// Annotation bit field
        /// </summary>
        IList<UInt16> Annotation { get; }

        /// <summary>
        /// Channel of the observed pulse. Channels are number starting at 0,
        /// ordered by ascending emission wavelength
        /// </summary>
        IList<byte> Channel { get; }

        /// <summary>
        /// Frame of movie when pulse began
        /// </summary>
        IList<float> StartFrame { get; }

        /// <summary>
        /// Duration of pulse in frames
        /// </summary>
        IList<float> WidthInFrames { get; }

        /// <summary>
        /// The maximum dye signal, in normalized [0,1] units
        /// </summary>
        IList<float> MaxSignal { get; }

        /// <summary>
        /// Number of wraps around the template sequence
        /// </summary>
        IList<UInt16> NumWraps { get; }

        /// <summary>
        /// Position on the template sequence
        /// </summary>
        IList<UInt16> Position { get; }
    }

    /// <summary>
    /// The base ZMW events interface
    /// </summary>
    public interface IZmwEvents
    {
        // Can promote this, in lieu of "NumBases", "NumPulses", etc.
        // but a lot of files will be hit with changes.
        //
        //uint Count { get; }
        
        /// <summary>
        /// A reference to the ISequencingZmw object that give basic metadata about the Zmw
        /// </summary>
        ISequencingZmw Zmw { get; }
    }

    /// <summary>
    /// The basecalls from a single ZMW. PulseIndex is an index into the pulse feature arrays of the corresponding
    /// IZwmPulses.  This information will only be available after a base-calling phase.
    /// </summary>
    public interface IZmwBases : IZmwEvents
    {
        /// <summary>
        /// Number of observed base calls
        /// </summary>
        uint NumBases { get; }

        /// <summary>
        /// Base identity observed.  ASCII code representing the observed base.
        /// </summary>
        IList<char> Base { get; }
        
        /// <summary>
        /// Phred-like quality value of the base
        /// </summary>
        IList<byte> QV { get; }
        
        /// <summary>
        /// An index into the pulse stream corresponding to the pulse that was called as this base
        /// </summary>
        IList<int> PulseIndex { get; }
        
        /// <summary>
        /// The called base sequence
        /// </summary>
        string Sequence { get; }

        #region Rich Quality Value Fields

        /// <summary>
        /// Boolean flag to tell if we have rich QV data or not
        /// </summary>
        bool HaveRichQvData { get; }

        /// <summary>
        /// Phred-style quality value indicating the probability that the current base is an insertion
        /// </summary>
        IList<byte> InsertionQV { get; }

        /// <summary>
        /// Phred-style quality value indicating the total probability of a deleted base before the current base,
        /// excluding the merge case.
        /// </summary>
        IList<byte> DeletionQV { get; }

        /// <summary>
        /// Phred-style quality value indicating the total probability that a merged pulse-call produced this base
        /// </summary>
        IList<byte> MergeQV { get; }

        /// <summary>
        /// The ASCII code of the most likely base to have been deleted before the current base.
        /// The ASCII code for 'N' represents a 'roughly constant' deletion probability over all bases. 
        /// </summary>
        IList<char> DeletionTag { get; }

        /// <summary>
        /// Phred-style quality value indicating the total probability that the current basecall is a substitution error
        /// </summary>
        IList<byte> SubstitutionQV { get; }

        /// <summary>
        /// The ASCII code of the most likely alternative basecall at this position
        /// </summary>
        IList<char> SubstitutionTag { get; }

        bool HaveKineticData { get; }

        IList<ushort> WidthInFrames { get; }
        IList<ushort> PreBaseFrames { get; }

        #endregion

        /// <summary>
        /// Metrics that are derived from these basecalls
        /// </summary>
        IZmwMetricsBases Metrics { get; }
    }

    /// <summary>
    /// An extension of IZmwBases to represent base calls that are the result of a Circular Consensus
    /// procedure. The PulseIndex and Metrics fields are no longer meaningful and should not be used.
    /// The extra properties defined here indicate the regions of the raw basecalls stream that were
    /// used in the ciruclar consensus
    /// </summary>
    public interface IZmwConsensusBases : IZmwBases
    {
        /// <summary>
        /// Regions of raw basecalls that were used during consensus calling
        /// </summary>
        DelimitedSeqReg[] InsertRegions { get; }

        /// <summary>
        /// Total number of passes used, including partial passes
        /// </summary>
        int NumPasses { get; }

        /// <summary>
        /// Predicted average accuracy of consensus sequence
        /// </summary>
        float PredictedAccuracy { get; }

        /// <summary>
        /// How many mutations distinguish this template from the final template?
        /// </summary>
        /// <value>The number of mutations.</value>
        int NumberOfMutations { get; }

        /// <summary>
        /// How many mutations were attempted?
        /// </summary>
        /// <value>The number of tried mutations.</value>
        int NumberOfTriedMutations { get; set; }


        int IterationsTaken { get; set; }

        /// <summary>
        /// How much processing time did the consensus generation take in seconds
        /// </summary>
        /// <value>The processing time.</value>
        double ProcessingTime { get; set; }

        ///// <summary>
        ///// Length of the insert read
        ///// </summary>
        //int InsertReadLength { get; }
    }

    /// <summary>
    /// Outcome specification of the original "Productivity" metric, as the (symantically problematic) 0/1/2.
    /// </summary>
    public enum ProductivityClass
    {
        Empty,
        Productive,
        Other,
        
        NotDefined = 255
        // Storage limited
    };

    /// <summary>
    /// Outcome specification of a finer-grained productivity metric; see bug 19274.
    /// </summary>
    public enum LoadingClass
    {
        Empty,
        FullHqRead,
        PartialHqRead,
        Multiload,
        Indeterminate,

        NotDefined = 255
        // Storage limited
    };

    /// <summary>
    /// Base-related summary metrics.
    /// These metrics are "storage optional" -- providers that are backed by persistent storage
    /// should handle these properties with methods that fall back to on-the-fly computation
    /// when the data are not present.
    /// </summary>
    public interface IZmwMetricsBases
    {
        IList<RegionAnnotator.Region> Regions { get; }

        /// <summary>
        /// The fraction of bases called by channel
        /// </summary>
        float[] BaseFraction { get; }

        /// <summary>
        /// The average (global) pulse rate of called bases in the ZMW (1/sec)
        /// </summary>
        float BaseRate { get; }

        /// <summary>
        /// The mean pulse width of called bases in the ZMW (sec)
        /// </summary>
        float BaseWidth { get; }

        /// <summary>
        /// The robust mean pulse IPD of called bases in the ZMW (sec)
        /// </summary>
        float BaseIpd { get; }

        /// <summary>
        /// The percentage of IPDs in HQRegions above a threshold
        /// </summary>
        float Pausiness { get; }

        /// <summary>
        /// An estimate of the local base rate, i.e. the polymerase
        /// rate excluding pause events.
        /// </summary>
        float LocalBaseRate { get; }

        /// <summary>
        /// An estimate of the dark base rate, i.e. the polymerase
        /// rate excluding pause events, that is expected when the chip isn't illuminated
        /// </summary>
        float DarkBaseRate { get; }

        /// <summary>
        /// A score corresponding to predicted accuracy of the read
        /// </summary>
        float ReadScore { get; }

        /// <summary>
        /// A classification of ZMW productivity
        /// </summary>
        ProductivityClass Productivity { get; }

        /// <summary>
        /// A classification of ZMW loading state
        /// </summary>
        LoadingClass Loading { get; }

        /// <summary>
        /// Control read accuracy 
        /// </summary>
        float ControlReadQual { get; } 

        /// <summary>
        /// Control read length 
        /// </summary>
        uint ControlReadLen { get; } 

        /// <summary>
        /// The SNR of the pulses inside the HQRegion
        /// </summary>
        float[] HQRegionSNR { get;  }

        /// <summary>
        /// The intra-pulse signal std of pulses inside the HQRegion
        /// </summary>
        float[] HQRegionIntraPulseStd { get; }

        /// <summary>
        /// The SNR of the pulses inside the HQRegion
        /// </summary>
        float HQRegionStartTime { get; }

        /// <summary>
        /// The SNR of the pulses inside the HQRegion
        /// </summary>
        float HQRegionEndTime { get; }

        /// <summary>
        /// Channel-mean base quality value
        /// Mean phred-style QV by base channel over the read
        /// </summary>
        float[] CmBasQv { get; }

        /// <summary>
        /// Channel-mean insertion quality value
        /// Mean InsertionQV by base channel over the read
        /// </summary>
        float[] CmInsQv { get; }

        /// <summary>
        /// Channel-mean deletion quality value
        /// Mean DeletionQV by base channel over the read
        /// </summary>
        float[] CmDelQv { get; }

        /// <summary>
        /// Channel-mean substitution quality value
        /// Mean SubstitutionQV by base channel over the read
        /// </summary>
        float[] CmSubQv { get; }

        /// <summary>
        /// Read-mean base quality value
        /// Mean phred-style QV over all bases in the read
        /// </summary>
        float RmBasQv { get; }

        /// <summary>
        /// Read-mean insertion quality value
        /// Mean InsertionQV over all bases in the read
        /// </summary>
        float RmInsQv { get; }

        /// <summary>
        /// Read-mean deletion quality value
        /// Mean DeletionQV over all bases in the read
        /// </summary>
        float RmDelQv { get; }

        /// <summary>
        /// Read-mean substitution quality value
        /// Mean SubstitutionQV over all bases in the read
        /// </summary>
        float RmSubQv { get; }

        /// <summary>
        /// Mean estimate of pkmid (SignalModel intra-pulse amplitude) in HQ region
        /// </summary>
        float[] HQRegionEstPkmid { get; }

        /// <summary>
        /// Mean estimate of pkstd (SignalModel intra-pulse sigma) in HQ region
        /// </summary>
        float[] HQRegionEstPkstd { get; }

        /// <summary>
        /// (Observed intra-pulse variance) / (Model-predicted intra-pulse variance) in the HQ region.
        /// The model prediction uses the observed mean pkmid in the HQ region as the pulse amplitude.
        /// </summary>
        float[] HQRegionPkzvar { get; }

        #region Time-Dependent Metrics

        /// <summary>
        /// The interval in sec over which all time-dependent metrics are computed.
        /// </summary>
        float TdmIntervalSec { get; }

        /// <summary>
        /// The cutoff for designating an HQ base IPD as a pause event.
        /// </summary>
        float PauseIpdThreshSec { get; }

        /// <summary>
        /// The number of base IPD values that are above the pause threshold,
        /// in HQ regions, and by fixed time interval along the read.  This
        /// contributes to the numerator for a PauseFraction calculation.
        /// </summary>
        short[] NumPauseVsT { get; }

        /// <summary>
        /// The total number of bases in HQ regions, by fixed time intervals
        /// along the read.  This contributes to the denominator for a
        /// PauseFraction calculation.  
        /// </summary>
        short[] NumBaseVsT { get; }

        /// <summary>
        /// The base incorporation rate, by fixed time intervals along the read.
        /// </summary>
        float[] BaseRateVsT { get; }

        #endregion
    }

    public enum TraceRepresentation
    {
        Camera,
        DyeWeightedSum,
        Multicomponent,
    }


    /// <summary>
    /// A class to hold a variety of Extension methods that 
    /// make it easier to work with Primary Analysis data.
    /// </summary>
    public static class AnalysisDataHelpers
    {
        public static RegionAnnotator.Region HQRegion(this IZmwBases b)
        {
            return b.Metrics.Regions.Where(r => r.Type.Type == "HQRegion").FirstOrDefault();
        }

        public static IEnumerable<RegionAnnotator.Region> HQRegions(this IZmwBases b)
        {
            return b.Metrics.Regions.Where(r => r.Type.Type == "HQRegion");
        }

        public static IEnumerable<RegionAnnotator.Region> AdapterHits (this IZmwBases b)
        {
            return b.Metrics.Regions.Where(r => r.Type.Type == "Adapter");
        }


        public static int MeanInsert(this IZmwBases b)
        {
            var hits = b.AdapterHits().ToArray();

            if(hits.Length >= 2)
            {
                var lSum = 0;

                var s0 = hits[0].End;

                for (int i = 1; i < hits.Length; i++)
                {
                    lSum += hits[i].Start - s0;
                    s0 = hits[i].End;
                }

                return lSum/(hits.Length - 1);
            }

            return -1;
        }

        public static int RegionLength(this IZmwBases b)
        {
            var region = b.HQRegion();

            if(region != null && b.Metrics.Productivity == ProductivityClass.Productive)
            {
                return region.Length;
            }

            return (int) b.NumBases;
        }

        public static float[] BaseStartTime(this IZmwBases b)
        {
            var fr = b.Zmw.Movie.FrameRate;
            var bw = b.WidthInFrames.Select(v => v / fr).ToArray();
            var ipd = b.PreBaseFrames.Select(v => v / fr).ToArray();

            var t1 = new float[bw.Length];
            float lastTime = 0;

            for (int i = 0; i < bw.Length; i++)
            {
                t1[i] = lastTime + ipd[i];
                lastTime = t1[i] + bw[i];
            }

            return t1;
        }
    }    
}
