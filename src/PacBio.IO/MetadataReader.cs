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
using System.Text.RegularExpressions;
using System.Xml;
using System.Xml.Linq;
using System.Xml.XPath;
using PacBio.IO.ICD;
using PacBio.HDF;
using PacBio.Utils;

// ReSharper disable EmptyGeneralCatchClause
namespace PacBio.IO
{
    [HdfIcdEntry(Path = "/ScanData/DyeSet/Analog[0]", Detail = "List of analogs, i=0…nA, as individual groups")]
    class HDFAnalogSpec : IAnalogSpec
    {
        public class Icd : HdfIcd<HDFAnalogSpec>
        {
            public Icd(bool flatten = false) : base(Mine(), flatten) { }
            private static IEnumerable<HdfIcdEntry> Mine()
            {
                return new[]
                {
                    // Explicit entries
                    new HdfIcdEntry{Path = "Type", Detail = "Name of SpectralCal record type", Deprecated = true},
                    new HdfIcdEntry{Path = "TypeId", Detail = "Id of SpectralCal record type", Deprecated = true},
                    new HdfIcdEntry{Path = "RecordName", Detail = "Name of the record that produced the calibration data", Deprecated = true},
                    new HdfIcdEntry{Path = "RecordDate", Detail = "Date associated with RecordName", Deprecated = true},
                    new HdfIcdEntry{Path = "CalMovieDate", Detail = "Date of the spectral calibration acquisition", Deprecated = true},
                };
            }
        }

        public HDFAnalogSpec(IGroup analogSpec)
        {
            analogSpecDataset = analogSpec;

            bse = analogSpec.GetAttribute("Base").ReadSingleton<string>()[0];
            label = analogSpec.GetAttribute("Label").ReadSingleton<string>();
            wavelength = analogSpec.GetAttribute("Wavelength").ReadSingleton<float>();
            nucleotide = analogSpec.GetAttribute("Nucleotide").ReadSingleton<string>();
        }

        #pragma warning disable 0414
        private IGroup analogSpecDataset;
        #pragma warning restore 0414

        private char bse = 'N';
        [HdfIcdEntry(Detail = "DNA base associated with this analog: A, C, G, or T")]
        public char  Base
        {
	        get { return bse; }
        }

        private float wavelength;
        [HdfIcdEntry(Units = "nm", Detail = "Characteristic dye emission wavelength")]
        public float  Wavelength
        {
	        get { return wavelength; }
        }

        private string label = "";
        [HdfIcdEntry(Detail = "Name of the dye label")]
        public string  Label
        {
	        get { return label; }
        }

        private string nucleotide = "";
        [HdfIcdEntry(Detail = "Nucleotide [and linker] detailed identifier", Deprecated = true)]
        public string  Nucleotide
        {
	        get { return nucleotide; }
        }

        private double concentration = Double.NaN;
        public double  Concentration
        {
	        get { return concentration; }
        }
    }

    /// <summary>
    /// A factory class to return the correct type of MetadataReader as IZmwSource.
    /// </summary>
    public static class MetadataReader
    {
        public static IZmwSource GetMetadataReader(IGroup scanData, IGroup zmwGroup)
        {
            // Key the metadata reader on the PlatformID field in ScanData/RunInfo/PlatformId
            // If that field doesn't exist we default to Astro
            var runInfoGroup = (IGroup)(scanData.GetChild("RunInfo"));
            var platformAttr = runInfoGroup.GetAttribute("PlatformId");

            if (platformAttr == null)
                return new AstroMetadataReader(scanData, zmwGroup);

            var platformId = platformAttr.ReadSingleton<uint>();

            if (platformId == 1)
                return new AstroMetadataReader(scanData, zmwGroup);

            if (platformId == 2)
                return new SpringfieldMetadataReader(scanData, zmwGroup);

            if (platformId == 3)
                return new SequelPocMetadataReader(scanData, zmwGroup);

            throw new Exception("Unrecognized PlatformId in ScanData/RunInfo");
        }
    }

    /// <summary>
    /// Contains one entry for every ZMW in the physical grid. 
    /// ByZmwIndex is the same as direct indexing
    /// </summary>
    [HdfIcdEntry(Path = "/ScanData", Detail = "Movie acquisition metadata")]
    public class SpringfieldMetadataReader : DataSource<ISequencingZmw>, IZmwSource
    {
        #region Members

        private readonly IGroup scanDataGroup;
        private readonly string filename;

        private string movieName = "";
        private string softwareVersion;
        private string changelist;
        private string bindingKit;
        private string sequencingKit;
        private string baseCallerChangelistID;
        private DateTime dateCreated;

        #pragma warning disable 0414
        private IAnalogSpec[] analogs;
        #pragma warning restore 0414

        private ISequencingZmw[] sequencingZmws;
        private readonly int numZmws;

        #endregion

        #region ICD

        public class Icd : HdfIcd<SpringfieldMetadataReader>
        {
            public Icd(bool flatten = false) : base(Mine(), flatten)
            {
                Concat(new HDFAnalogSpec.Icd(flatten));
            }
            private static IEnumerable<HdfIcdEntry> Mine()
            {
                return new[]
                    {
                        // Explicit entries
                        new HdfIcdEntry {Path = "FormatVersion", Detail = "ID/version of the ScanData schema"},
                        new HdfIcdEntry {Path = "SoftwareVersion", Detail = "ID/version of the software that created the ScanData group"},
                        new HdfIcdEntry {Path = "ChangeListID", Detail = "Revision ID of the software that created the ScanData group"},
                        new HdfIcdEntry {Path = "DateCreated", Units = "ISO 8601", Detail = "Time-stamp for creation of the ScanData group"},

                        // Acquisition params, including here metadata not supported by this API
                        new HdfIcdEntry {Path = "AcqParams", Detail = "Instrument acquisition parameters"},
                        new HdfIcdEntry {Path = "AcqParams/Look", Detail = "Set number of the acquisition"},
                        new HdfIcdEntry {Path = "AcqParams/LaserOnFrameValid", Detail = "Non-zero if LaserOn frame estimate was available"},
                        new HdfIcdEntry {Path = "AcqParams/HotStartFrameValid", Detail = "Non-zero if HotStart frame estimate was available"},
                        new HdfIcdEntry {Path = "AcqParams/NumLasers", Detail = "The number of excitation sources, nL"},
                        new HdfIcdEntry {Path = "AcqParams/NumCameras", Detail = "The number of acquisition cameras"},
                        new HdfIcdEntry {Path = "AcqParams/WhenStarted", Units = "ISO 8601", Detail = "Time-stamp for start of acquisition"},

                        // Chip metadata...
                        new HdfIcdEntry {Path = "ChipInfo", Detail = "Chip type and layout metadata"},
                        new HdfIcdEntry {Path = "ChipInfo/LayoutName", Detail = "Identifier or name of the chip layout"},
                        new HdfIcdEntry {Path = "ChipInfo/LayoutId", Detail = "Identifier of the chip layout", Deprecated = true},
                        new HdfIcdEntry {Path = "ChipInfo/NumLooks", Detail = "Number of looks (sets) on the chip layout"},
                        new HdfIcdEntry {Path = "ChipInfo/GridDims", Detail = "X/Y grid dimensions of the chip layout"},
                        new HdfIcdEntry {Path = "ChipInfo/ChipId", Detail = "Unique identifier of this chip, as <strip_serial_number>-<index>", Deprecated = true},
                        new HdfIcdEntry {Path = "ChipInfo/NumDiffusion", Detail = "The number of ZMWs available for diffusion-background estimation"},
                        new HdfIcdEntry {Path = "ChipInfo/DiffusionXY", Detail = "Chip X/Y coordinates of the ZMWs available for diffusion-background estimation"},

                        // Dye-set metadata
                        new HdfIcdEntry {Path = "DyeSet", Detail = "Dye-labeled analog metadata"},
                        new HdfIcdEntry {Path = "DyeSet/Name", Detail = "Dye-set name", Deprecated = true},
                        new HdfIcdEntry {Path = "DyeSet/BaseMap", Detail = "Mapping of dyes (ordered by wavelength) to bases"},
                        new HdfIcdEntry {Path = "DyeSet/NumAnalog", Detail = "Number of analogs (dyes) in the dyeset"},

                        // Chemistry protocol (deprecated)
                        new HdfIcdEntry {Path = "Protocol", Detail = "Chemistry protocol metadata"},
                        new HdfIcdEntry {Path = "Protocol/Name", Detail = "Name of run protocol", Deprecated = true},
                        new HdfIcdEntry {Path = "Protocol/InsertSize", Detail = "Approximate size of DNA insert", Deprecated = true},

                        // General run inforamtion
                        new HdfIcdEntry {Path = "RunInfo", Detail = "General run, movie and configuration data"},
                        new HdfIcdEntry {Path = "RunInfo/PlatformId", Detail = "Platform identifier", Units = "1:Astro,2:PacBio-RS,3:SequelAstro"},
                        new HdfIcdEntry {Path = "RunInfo/InstrumentName", Detail = "Instrument instance identifier"},
                        new HdfIcdEntry {Path = "RunInfo/MovieName", Detail = "Movie context string"},
                        new HdfIcdEntry {Path = "RunInfo/RunCode", Detail = "?? Definition ??"},
                        new HdfIcdEntry {Path = "RunInfo/BindingKit", Detail = "Binding kit version"},
                        new HdfIcdEntry {Path = "RunInfo/SequencingKit", Detail = "Sequencing kit version"},
                        new HdfIcdEntry {Path = "RunInfo/PlatformName", Detail = "Platform name", Deprecated = true},
                        new HdfIcdEntry {Path = "RunInfo/InstrumentId", Detail = "Unique instrument ID", Deprecated = true},
                        new HdfIcdEntry {Path = "RunInfo/RunId", Detail = "Unique acquisition ID", Deprecated = true}
                    };
            }
        }

        #endregion

        #region Properties

        public uint NumChannels { get; private set; }

        private readonly ZmwIndexer indexer;
        public ZmwIndexer ZmwIndexer { get { return indexer; } }

        // Provided to facilitate group copy between file types.
        public IGroup ScanDataGroup
        {
            get { return scanDataGroup; }
        }

        // This is its own ZmwSource, backing the DataSource<ISequencingZmw> implementation.
        public override IZmwSource ZmwSource
        {
            get
            {
                return this;
            }
            protected set
            {
                // Don't use the set in this case.
                throw new NotImplementedException();
            }
        }

        #endregion

        #region Structors

        public SpringfieldMetadataReader(IGroup scanData, IGroup zmwGroup)
        {
            filename = scanData.File.Name;

            if (zmwGroup.GetChild("HoleXY") == null)
                throw new Exception("No a valid HDF5 ZMW group");
            
            indexer = new ZmwIndexer(zmwGroup);
            numZmws = indexer.NumZmws;

            scanDataGroup = scanData;
            Load();
        }

        private void Load()
        {
            LoadAcqParams();
            LoadChipInfo();
            LoadRunInfo();

            IGroup dyeSet = scanDataGroup.CreateGroup("DyeSet");

            var verAttr = scanDataGroup.GetAttribute("SoftwareVersion") ?? scanDataGroup.GetAttribute("ProviderName");
            softwareVersion = verAttr != null ? verAttr.ReadSingleton<string>() : "Unknown version";
            
            var clAttr = scanDataGroup.GetAttribute("ChangeListID");
            changelist = clAttr != null ? clAttr.ReadSingleton<string>() : "Unknown version";

            dateCreated = new DateTime(1900, 1, 1);

            try
            {
                var dateAttribute = scanDataGroup.GetAttribute("DateCreated");
                
                if(dateAttribute != null)
                {
                    var dateString = dateAttribute.ReadSingleton<string>();
                    if (!DateTime.TryParse(dateString, out dateCreated))
                        dateCreated = new DateTime(1900, 1, 1);
                }
            }            
            catch (Exception)           
            {
            }

            // Make an IAnalogSpec for each of the analogs present
            NumChannels = dyeSet.GetAttribute("NumAnalog").ReadSingleton<UInt16>();
            analogs = NumChannels.Fill<IAnalogSpec>(i => new HDFAnalogSpec(dyeSet.CreateGroup(String.Format("Analog[{0}]", i))));

            // Setup our actual sequencing ZMWs
            sequencingZmws = numZmws.Fill(i => new HDFSequencingZMW(this, i));
        }

        #endregion

        #region AcqParam Data

        //private byte look;
        private float frameRate;
        private uint numFrames;
        private int laserOnFrame;
        private int hotStartFrame;

        private void LoadAcqParams()
        {
            var acqPars = (IGroup) scanDataGroup.GetChild("AcqParams");

            // Make sure we have acqParams
            if (acqPars == null)
                throw new Exception("Bad Trace file");

            // Get basic movie data
            //look = acqPars.GetAttribute("Look").ReadSingleton<byte>();
            frameRate = acqPars.GetAttribute("FrameRate").ReadSingleton<float>();
            numFrames = acqPars.GetAttribute("NumFrames").ReadSingleton<uint>();
            laserOnFrame = acqPars.GetAttribute("LaserOnFrame").ReadSingleton<int>();
            hotStartFrame = acqPars.GetAttribute("HotStartFrame").ReadSingleton<int>();

            var dyeSetGroup = (IGroup) scanDataGroup.GetChild("DyeSet");
            var bm = dyeSetGroup.GetAttribute("BaseMap").ReadSingleton<string>();
            BaseMap = bm.ToCharArray();
        }

        #endregion

        #region ChipInfo Data

        /*
        private string chipLayoutName;
        private uint chipLayoutId;
        private int[] gridDims;
        private byte numLooks;
        private string chipId;
        */

        private void LoadChipInfo()
        {
            /*
            IGroup chipInfo = scanDataGroup.CreateGroup("ChipInfo");
            chipLayoutName = chipInfo.GetAttribute("LayoutName").ReadSingleton<string>();
            chipLayoutId = chipInfo.GetAttribute("LayoutId").ReadSingleton<uint>();
            gridDims = (int[]) chipInfo.GetAttribute("GridDims").Read();
            numLooks = chipInfo.GetAttribute("NumLooks").ReadSingleton<byte>();
            chipId = chipInfo.GetAttribute("ChipId").ReadSingleton<string>();
            */ 
        }

        #endregion

        #region RunInfo Data

        public string PlatformName { get; private set; }
        public uint PlatformId { get; private set; }
        public string InstrumentName { get; private set; }
        public uint InstrumentId { get; private set; }
        public string RunCode { get; private set; }
        public uint RunId { get; private set; }
        public char[] BaseMap { get; private set; }

        private void LoadRunInfo()
        {
            var runInfo = scanDataGroup.CreateGroup("RunInfo");

            PlatformName = runInfo.GetAttribute("PlatformName").ReadSingleton<string>();
            PlatformId = runInfo.GetAttribute("PlatformId").ReadSingleton<uint>();
            InstrumentName = runInfo.GetAttribute("InstrumentName").ReadSingleton<string>();
            InstrumentId = runInfo.GetAttribute("InstrumentId").ReadSingleton<uint>();

            RunCode = runInfo.GetAttribute("RunCode").ReadSingleton<string>();
            RunId = runInfo.GetAttribute("RunId").ReadSingleton<uint>();

            // load the sequencing chemistry information, if it's available from the hdf5,
            // otherwise load it from the metadata.xml
            var bkAttr = runInfo.GetAttribute("BindingKit");
            bindingKit = bkAttr != null ? bkAttr.ReadSingleton<string>() : null;

            var skAttr = runInfo.GetAttribute("SequencingKit");
            sequencingKit = skAttr != null ? skAttr.ReadSingleton<string>() : null;

            var baseCalls = scanDataGroup.File.CreateGroup("PulseData").CreateGroup("BaseCalls");
            var bcAttr = baseCalls.GetAttribute("ChangeListID");
            baseCallerChangelistID = bcAttr != null ? bcAttr.ReadSingleton<string>() : null;

            try
            {
                movieName = runInfo.GetAttribute("MovieName").ReadSingleton<string>();              
            } 
            catch(Exception)
            {                
            }
        }

        private string XmlSelect(XPathNavigator nav, XmlNamespaceManager man, string xpath)
        {
            return nav.SelectSingleNode(xpath, man).InnerXml;
        }

        private void LoadMetadataXml()
        {
            // if we're missing even one value, load everything from the metadata xml
            if (bindingKit == null || sequencingKit == null || baseCallerChangelistID == null)
            {
                var up = Path.GetDirectoryName(Path.GetDirectoryName(Path.GetFullPath(FileName)));
                var metadataLocation = Path.Combine(up, MovieName + ".metadata.xml");

                if (!System.IO.File.Exists(metadataLocation))
                    throw new PacBio.Data.ChemistryLookupException(String.Format("Could not find, {0}", metadataLocation));

                var doc = new XPathDocument(metadataLocation);
                var nav = doc.CreateNavigator();
                var man = new XmlNamespaceManager(nav.NameTable);

                man.AddNamespace("pb", "http://pacificbiosciences.com/PAP/Metadata.xsd");

                try
                {
                    bindingKit = XmlSelect(nav, man, "pb:Metadata/pb:BindingKit/pb:PartNumber");
                    sequencingKit = XmlSelect(nav, man, "pb:Metadata/pb:SequencingKit/pb:PartNumber");
                    baseCallerChangelistID = XmlSelect(nav, man, "pb:Metadata/pb:InstCtrlVer");
                }
                catch
                {
                    throw new PacBio.Data.ChemistryLookupException(
                        String.Format("Could not extract chemistry information from '{0}'", metadataLocation));
                }
            }
        }

        #endregion

        #region IList
  
        public override ISequencingZmw this[int index]
        {
            get
            {
                return sequencingZmws[index];
            }
            set
            {
                throw new NotImplementedException();
            }
        }

        #endregion

        #region ICollection

        public override bool Contains(ISequencingZmw item)
        {
            return sequencingZmws.Contains(item);
        }

        public override int Count
        {
            get { return sequencingZmws.Length; }
        }

        #endregion

        #region IEnumerable

        public override IEnumerator<ISequencingZmw> GetEnumerator()
        {
            return sequencingZmws.AsEnumerable().GetEnumerator();
        }

        #endregion

        #region IDisposable

        protected override void Dispose(bool disposing)
        {}

        #endregion

        #region ISourceIdentifier

        public override string FileName
        {
            get { return filename; }
        }

        public override string SoftwareVersion
        {
            get { return softwareVersion; }
        }

        public override string ChangelistID
        {
            get { return changelist; }
        }

        public override DateTime DateCreated
        {
            get { return dateCreated; }
        }

        public string BindingKit
        {
            get
            {
                LoadMetadataXml();

                return bindingKit;
            }
        }

        public string SequencingKit
        {
            get
            {
                LoadMetadataXml();

                return sequencingKit;
            }
        }

        public string BaseCallerChangelistID
        {
            get
            {
                LoadMetadataXml();

                return baseCallerChangelistID;
            }
        }

        /// <summary>
        /// Return the sequencing enzyme/chemistry/software-version as a 3-tuple,
        /// usually for the purposes of determining the sequencing chemistry (e.g P4-C2).
        /// Raises ChemistryLookupError in the event of an error or missing metadata.xml
        /// </summary>
        /// <value>The chemistry barcode triple.</value>
        public PacBio.Data.ChemistryBarcodeTriple ChemistryBarcode
        {
            get
            {
                return new PacBio.Data.ChemistryBarcodeTriple(BindingKit, SequencingKit, BaseCallerChangelistID);
            }
        }

        /// <summary>
        /// Return the Quiver Sequencing Chemistry string (e.g. P4-C2)
        /// </summary>
        /// <value>The sequencing chemistry.</value>
        public string SequencingChemistry
        {
            get
            {
                return "P6-C4";
                return PacBio.Data.Chemistry.DecodeTriple(ChemistryBarcode);
            }
        }

        #endregion

        #region IDataSource

        /// <summary>
        /// Overrides the base-class implementation for efficiency 
        /// </summary>
        /// <param name="holeNum"></param>
        /// <returns></returns>
        public override ISequencingZmw ByHoleNumber(int holeNum)
        {
            return sequencingZmws[GetIndexByHoleNumber(holeNum)];
        }

        public override ISequencingZmw ByXY(int x, int y)
        {
            try
            {
                return sequencingZmws[GetIndexByHoleXY(x, y)];
            }
            catch (KeyNotFoundException)
            {
                throw new Exception(String.Format("ZMW X:{0}, Y{1} does not exist in this file", x, y));
            }
        }

        #endregion

        private class HDFSequencingZMW : ISequencingZmw
        {
            private readonly SpringfieldMetadataReader reader;

            public HDFSequencingZMW(SpringfieldMetadataReader reader, int index)
            {
                this.reader = reader;
                this.index = index;
            }

            private readonly int index;

            public IMovieMetadata Movie
            {
                get { return reader; }
            }

            public int HoleNumber
            {
                get { return reader.indexer.HoleNumber[index]; }
            }

            // Astro defaults -- get more fancy later
            public int X
            {
                get { return reader.indexer.HoleXY[index].X; }
            }

            public int Y
            {
                get { return reader.indexer.HoleXY[index].Y; }
            }

            public Int16 HoleChipLook
            {
                get { return reader.indexer.HoleChipLook[index];  }
            }

            // Handle the Astro cases - need to Expand in the future.
            public ZmwType ZmwType
            {
                get { return reader.indexer.ZmwType[index]; }
            }

            public override string ToString()
            {
                return String.Format("{0}/{1}", Movie.MovieName, HoleNumber);
            }
        }



        #region IMovieMetadata

        public string MovieName
        {
            get { return movieName; }
        }

        [HdfIcdEntry(Path = "AcqParams/", Units = "Hz", Detail = "Acquisition frame rate")]
        public float FrameRate
        {
            get { return frameRate; }
        }

        [HdfIcdEntry(Path = "AcqParams/", Detail = "Number of frames acquired")]
        public uint NumFrames
        {
            get { return numFrames; }
        }

        [HdfIcdEntry(Path = "AcqParams/", Detail = "Frame number estimate of a laser-on event; may be <= 0")]
        public int LaserOnFrame
        {
            get { return laserOnFrame; }
        }

        [HdfIcdEntry(Path = "AcqParams/", Detail = "Frame number estimate of a hot start event; may be <= 0")]
        public int HotStartFrame
        {
            get { return hotStartFrame; }
        }

        #endregion

        #region IZmwIndexer Members

        public int GetIndexByHoleNumber(int holeNumber)
        {
            return indexer.GetIndexByHoleNumber(holeNumber);
        }

        public int GetIndexByHoleXY(int x, int y)
        {
            return indexer.GetIndexByHoleXY(x, y);
        }
        
        #endregion
    }
    // ReSharper restore EmptyGeneralCatchClause


    /// <summary>
    /// Contains one entry for every ZMW in the physical grid. 
    /// ByZmwIndex is the same as direct indexing
    /// </summary>
    public class AstroMetadataReader : DataSource<ISequencingZmw>, IZmwSource
    {
        #region Members

        private IGroup scanDataGroup;
        private float frameRate;
        private uint numFrames;
        private int laserOnFrame;
        private int hotStartFrame;
        private float[] laserIntensity;
        private float[] laserPower;
        private float cameraBias;
        private float cameraBiasStd;
        private float cameraGain;
        private float aduGain;

        private string runCode = "";
        private string movieName = "";
        private string instrumentName = "";
        private string softwareVersion;
        private string changelist;
        private DateTime dateCreated;

        private IAnalogSpec[] analogs;
        private sbyte[,] chipMask;
        private ISequencingZmw[] sequencingZmws;

        #endregion

        #region Properties

        public uint NumChannels { get; private set; }
        public uint NumLines { get; private set; }
        public uint HolesPerLine { get; private set; }

        private ZmwIndexer indexer;
        public ZmwIndexer ZmwIndexer { get { return indexer; } }

        // Provided to facilitate group copy between file types.
        public IGroup ScanDataGroup
        {
            get { return scanDataGroup; }
        }

        // This is its own ZmwSource, backing the DataSource<ISequencingZmw> implementation.
        public override IZmwSource ZmwSource
        {
            get
            {
                return this;
            }
            protected set
            {
                // Don't use the set in this case.
                throw new NotImplementedException();
            }
        }

        #endregion

        #region Structors

        public AstroMetadataReader(string hdfFile)
        {
            // Open the file
            var chunks = HDFFile.Open(hdfFile, FileMode.Open, FileAccess.Read);

            // Load from the ScanData group
            var scanData = chunks.CreateGroup("ScanData");
            Load(scanData, null);
        }

        public AstroMetadataReader(IGroup scanData, IGroup zmwGroup)
        {
            Load(scanData, zmwGroup);
        }

        public void Load(IGroup scanData, IGroup zmwGroup)
        {
            scanDataGroup = scanData;

            var acqPars = (IGroup)scanDataGroup.GetChild("AcqParams");

            // Load AcqParams
            if (acqPars != null)
            {
                // Get basic movie data
                frameRate = acqPars.GetAttribute("FrameRate").ReadSingleton<float>();
                numFrames = acqPars.GetAttribute("NumFrames").ReadSingleton<uint>();
                laserOnFrame = acqPars.GetAttribute("LaserOnFrame").ReadSingleton<int>();
                hotStartFrame = acqPars.GetAttribute("HotStartFrame").ReadSingleton<int>();
                laserIntensity = (float[])acqPars.GetAttribute("LaserIntensity").Read();
                cameraGain = acqPars.GetAttribute("CameraGain").ReadSingleton<float>();
                aduGain = acqPars.GetAttribute("AduGain").ReadSingleton<float>();

                // These fields may have been added later.  This try/catch can be removed in late 2009.
                try
                {
                    laserPower = (float[])acqPars.GetAttribute("LaserPower").Read();
                    cameraBias = acqPars.GetAttribute("CameraBias").ReadSingleton<float>();
                    cameraBiasStd = acqPars.GetAttribute("CameraBiasStd").ReadSingleton<float>();
                }
                catch (Exception)
                {
                    laserPower = new[] { float.NaN, float.NaN };
                    cameraBias = -1.0f;
                    cameraBiasStd = -1.0f;
                }
            }
            else
            {
                frameRate = 100;
                numFrames = 0;
                laserOnFrame = 0;
                hotStartFrame = 0;
                laserIntensity = new[] { float.NaN, float.NaN };
                laserPower = new[] { float.NaN, float.NaN };
                cameraBias = -1.0f;
                cameraBiasStd = -1.0f;
                cameraGain = 500.0f;
                aduGain = 27.0f;
            }

            IGroup chipArray = scanDataGroup.CreateGroup("ChipArray");
            IGroup dyeSet = scanDataGroup.CreateGroup("DyeSet");
            IGroup runInfo = scanDataGroup.CreateGroup("RunInfo");

            // Get data out of the ChipArray group
            NumLines = chipArray.GetAttribute("NumLines").ReadSingleton<uint>();
            HolesPerLine = chipArray.GetAttribute("HolesPerLine").ReadSingleton<uint>();

            var chipMaskDataset = (IDataset)chipArray.GetChild("ChipMask");
            var cmo = chipMaskDataset.Read();
            chipMask = cmo as sbyte[,];

            // Get stuff out of the RunInfo group
            instrumentName = runInfo.GetAttribute("InstrumentName").ReadSingleton<string>();

            var verAttr = scanDataGroup.GetAttribute("SoftwareVersion") ?? scanDataGroup.GetAttribute("ProviderName");
            softwareVersion = verAttr != null ? verAttr.ReadSingleton<string>() : "Unknown version";

            var clAttr = scanDataGroup.GetAttribute("ChangeListID");
            changelist = clAttr != null ? clAttr.ReadSingleton<string>() : "Unknown version";

            dateCreated = Helpers.GetDateAttribute(scanDataGroup, "DateCreated");

            runCode = runInfo.GetAttribute("RunCode").ReadSingleton<string>();
            movieName = runInfo.GetAttribute("MovieName").ReadSingleton<string>();

            // Make an IAnalogSpec for each of the analogs present
            NumChannels = dyeSet.GetAttribute("NumAnalog").ReadSingleton<UInt16>();
            analogs = NumChannels.Fill<IAnalogSpec>(i => new HDFAnalogSpec(dyeSet.CreateGroup(String.Format("Analog[{0}]", i))));

            var bm = dyeSet.GetAttribute("BaseMap").ReadSingleton<string>();
            BaseMap = bm.ToCharArray();


            // Setup our actual sequencing ZMWs
            if (zmwGroup != null)
            {
                var holeNumber = ((uint[])((IDataset)zmwGroup.GetChild("HoleNumber")).Read()).Map(v => (int) v);
                sequencingZmws = holeNumber.Map(i => new AstroSequencingZmw(this, i));

                var xy = sequencingZmws.Map(s => new XYPoint { X = s.X, Y = s.Y });
                indexer = new ZmwIndexer(holeNumber, xy);
            }
            else
            {
                sequencingZmws = new ISequencingZmw[] { };
                indexer = new ZmwIndexer(new int[] { }, new XYPoint[] { });
            }
        }


        #endregion

        #region IList

        public override ISequencingZmw this[int index]
        {
            get
            {
                return sequencingZmws[index];
            }
            set
            {
                throw new NotImplementedException();
            }
        }

        #endregion

        #region ICollection

        public override bool Contains(ISequencingZmw item)
        {
            return sequencingZmws.Contains(item);
        }

        public override int Count
        {
            get { return sequencingZmws.Length; }
        }

        #endregion

        #region IEnumerable

        public override IEnumerator<ISequencingZmw> GetEnumerator()
        {
            return sequencingZmws.AsEnumerable().GetEnumerator();
        }

        #endregion

        #region IDisposable

        protected override void Dispose(bool disposing)
        {}

        #endregion

        #region ISourceIdentifier

        public override string SoftwareVersion
        {
            get { return softwareVersion; }
        }

        public override string ChangelistID
        {
            get { return changelist; }
        }

        public override DateTime DateCreated
        {
            get { return dateCreated; }
        }

        #endregion

        #region IDataSource

        public override ISequencingZmw ByHoleNumber(int holeNum)
        {
            return sequencingZmws[GetIndexByHoleNumber(holeNum)];
        }

        public override ISequencingZmw ByXY(int x, int y)
        {
            try
            {   // This applies to standard Astro layout only.        
                int holeNumber = (int)HolesPerLine * (x - 1) + (y - 1);
                return ByHoleNumber(holeNumber);
            }
            catch (Exception)
            {
                throw new Exception(String.Format("ZMW X:{0}, Y{1} does not exist in this file", x, y));
            }
        }

        #endregion

        class AstroSequencingZmw : ISequencingZmw
        {
            private readonly AstroMetadataReader reader;
            private readonly int holeNumber;

            public AstroSequencingZmw(AstroMetadataReader readerIn, int index)
            {
                reader = readerIn;
                holeNumber = index;
            }

            public IMovieMetadata Movie
            {
                get { return reader; }
            }

            public int HoleNumber
            {
                get { return holeNumber; }
            }

            // Astro defaults -- get more fancy later
            public int X
            {
                get { return (int)(holeNumber / reader.HolesPerLine) + 1; }
            }

            public int Y
            {
                get { return (int)(holeNumber % reader.HolesPerLine) + 1; }
            }

            public float XPlot
            {
                get { return 0; }
            }

            public float YPlot
            {
                get { return 0; }
            }

            public Int16 HoleChipLook
            {
                get { return 0; }
            }

            public float Radius
            {
                get
                {
                    return (float)Math.Sqrt(Math.Pow((15 * (float)(X - 16)), 2) + Math.Pow((5 * (float)(Y - 47)), 2));
                }
            }

            public float Theta
            {
                get { throw new NotImplementedException(); }
            }

            // Handle the Astro cases - need to Expand in the future.
            public ZmwType ZmwType
            {
                get
                {
                    if ((Y == 1 || Y > 90 || (Y >= 46 && Y <= 48 && X == 17)) || (X == 2 && Y == 2))
                        return ZmwType.InvalidData;

                    if (reader.chipMask[Y - 1, X - 1] == 1)
                        return ZmwType.AntiHole;

                    return ZmwType.Sequencing;
                }
            }
        }

        #region IMovieMetadata

        public string MovieName
        {
            get { return movieName; }
        }

        public string RunCode
        {
            get { return runCode; }
        }

        public string InstrumentName
        {
            get { return instrumentName; }
        }

        public uint InstrumentId
        {
            get { return (uint)(Math.Abs(instrumentName.GetHashCode()) % 1000); }
        }

        public float FrameRate
        {
            get { return frameRate; }
        }

        public uint NumFrames
        {
            get { return numFrames; }
        }

        public int LaserOnFrame
        {
            get { return laserOnFrame; }
        }

        public int HotStartFrame
        {
            get { return hotStartFrame; }
        }

        public float[] LaserIntensity
        {
            get { return laserIntensity; }
        }

        public float[] LaserPower
        {
            get { return laserPower; }
        }

        public sbyte[,] ChipMask
        {
            get { return chipMask; }
        }

        public IList<IAnalogSpec> Analogs
        {
            get { return analogs; }
        }

        public float CameraBias
        {
            get { return cameraBias; }
        }

        public float CameraBiasStd
        {
            get { return cameraBiasStd; }
        }

        public float CameraGain
        {
            get { return cameraGain; }
        }

        public float AduGain
        {
            get { return aduGain; }
        }

        public char[] BaseMap { get; private set; }
        
        public uint PlatformId
        {
            get { return 1; }
        }

        public string PlatformName
        {
            get { return "Astro"; }
        }

        public string BindingKit
        {
            get { return "Astro"; }
        }

        public string SequencingKit
        {
            get { return "Astro"; }
        }

        public string BaseCallerChangelistID
        {
            get { return "Astro"; }
        }

        public PacBio.Data.ChemistryBarcodeTriple ChemistryBarcode
        {
            get { throw new PacBio.Data.ChemistryLookupException("Astro chemistry information is unknown"); }
        }

        public string SequencingChemistry
        {
            get { return PacBio.Data.Chemistry.DecodeTriple(ChemistryBarcode); }
        }

        #endregion

        #region IZmwIndexer Members

        public int GetIndexByHoleNumber(int holeNumber)
        {
            return indexer.GetIndexByHoleNumber(holeNumber);
        }

        public int GetIndexByHoleXY(int x, int y)
        {
            return indexer.GetIndexByHoleXY(x, y);
        }

        #endregion
    }

}
