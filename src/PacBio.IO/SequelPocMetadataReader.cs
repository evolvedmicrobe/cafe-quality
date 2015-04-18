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
    /// <summary>
    /// Contains one entry for every ZMW in the physical grid. 
    /// ByZmwIndex is the same as direct indexing
    /// </summary>
    public class SequelPocMetadataReader : DataSource<ISequencingZmw>, IZmwSource
    {
        #region Members

        private IGroup scanDataGroup;
        private float frameRate;
        private uint numFrames;
        private int laserOnFrame;
        private int hotStartFrame;
        private float[] laserPower;
        private float cameraBias;
        private float cameraBiasStd;
        private float cameraGain;
        private float aduGain;

        private string movieName = "";
        private string instrumentName = "";
        private string softwareVersion;
        private string changelist;
        private DateTime dateCreated;

        private IAnalogSpec[] analogs;
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

        public SequelPocMetadataReader(string hdfFile)
        {
            // Open the file
            var chunks = HDFFile.Open(hdfFile, FileMode.Open, FileAccess.Read);

            // Load from the ScanData group
            var scanData = chunks.CreateGroup("ScanData");
            Load(scanData, null);
        }

        public SequelPocMetadataReader(IGroup scanData, IGroup zmwGroup)
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
                laserPower = new[] { float.NaN, float.NaN };
                cameraBias = -1.0f;
                cameraBiasStd = -1.0f;
                cameraGain = 500.0f;
                aduGain = 27.0f;
            }

            IGroup dyeSet = scanDataGroup.CreateGroup("DyeSet");
            IGroup runInfo = scanDataGroup.CreateGroup("RunInfo");

            // Get data out of the ChipArray group
            // SEQUEL PoC has no ChipArray/{NumLines,HolesPerLine}, ChipMask, LaserIntensity, or RunCode
            NumLines = 0;
            HolesPerLine = 1;

            // Get stuff out of the RunInfo group
            instrumentName = runInfo.GetAttribute("InstrumentName").ReadSingleton<string>();

            var verAttr = scanDataGroup.GetAttribute("SoftwareVersion") ?? scanDataGroup.GetAttribute("ProviderName");
            softwareVersion = verAttr != null ? verAttr.ReadSingleton<string>() : "Unknown version";

            var clAttr = scanDataGroup.GetAttribute("ChangeListID");
            changelist = clAttr != null ? clAttr.ReadSingleton<string>() : "Unknown version";

            dateCreated = Helpers.GetDateAttribute(scanDataGroup, "DateCreated");

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
                sequencingZmws = holeNumber.Map(i => new SequelPocSequencingZmw(this, i));

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
            {   // This applies to standard SequelPOC layout only.        
                int holeNumber = (int)HolesPerLine * (x - 1) + (y - 1);
                return ByHoleNumber(holeNumber);
            }
            catch (Exception)
            {
                throw new Exception(String.Format("ZMW X:{0}, Y{1} does not exist in this file", x, y));
            }
        }

        #endregion

        class SequelPocSequencingZmw : ISequencingZmw
        {
            private readonly SequelPocMetadataReader reader;
            private readonly int holeNumber;

            public SequelPocSequencingZmw(SequelPocMetadataReader readerIn, int index)
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

            // SequelPOC defaults -- get more fancy later
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

            // Handle the SequelPoC cases - need to Expand in the future.
            public ZmwType ZmwType
            {
                get { return ZmwType.Sequencing; }
            }
        }

        #region IMovieMetadata

        public string MovieName
        {
            get { return movieName; }
        }

        public string RunCode
        {
            get { return "SequelPoC"; }
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
            get { return new float[]{ 1.0f, 1.0f }; }
        }

        public float[] LaserPower
        {
            get { return laserPower; }
        }

        public sbyte[,] ChipMask
        {
            get { return new sbyte[,]{}; }
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
            get { return 3; }
        }

        public string PlatformName
        {
            get { return "SequelPoC"; }
        }

        public string BindingKit
        {
            get { return "SequelPoC"; }
        }

        public string SequencingKit
        {
            get { return "SequelPoC"; }
        }

        public string BaseCallerChangelistID
        {
            get { return "SequelPoC"; }
        }

        public PacBio.Data.ChemistryBarcodeTriple ChemistryBarcode
        {
            get { throw new PacBio.Data.ChemistryLookupException("SequelPOC chemistry information is unknown"); }
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
