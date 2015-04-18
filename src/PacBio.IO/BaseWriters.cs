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
using System.Reflection;
using PacBio.HDF;
using PacBio.IO.ICD;
using PacBio.Utils;

namespace PacBio.IO
{
    // Need to leave in some Type arguments to make Mono happy -- suppress reshaper warning
    // ReSharper disable RedundantTypeArgumentsOfMethod

    /// <summary>
    /// A temporary post-processing step -- we copy the PulseData\BasCcalls group from the .bas.h5 file to the .pls.h5 file
    /// for backward compatibility
    /// </summary>
    public static class BasecallsBackCopy
    {
        public static bool CopyHDFData(string destinationFile, string sourceFile, string subGroupName)
        {
            using (var targetFile = HDFFile.Open(destinationFile, FileMode.Open, FileAccess.ReadWrite))
            {
                using (var source = HDFFile.Open(sourceFile, FileMode.Open, FileAccess.Read))
                {
                    using (IGroup sourceBaseData = (IGroup)source.GetChild("PulseData"),
                                  destPulseData = (IGroup)targetFile.GetChild("PulseData"))
                    {
                        using (HDFIntraFile sourceItem = (HDFIntraFile)sourceBaseData.GetChild(subGroupName),
                                        targetItem = (HDFIntraFile)destPulseData.GetChild(subGroupName))
                        {
                            // Source item doesn't exist, can't copy
                            if (sourceItem == null)
                                return false;

                            // Clear out existing target item if it exists
                            if (targetItem != null)
                                destPulseData.RemoveChild(subGroupName);

                            sourceItem.Copy(destPulseData);
                            return true;
                        }
                    }
                }
            }
        }
    }


    public class WriterEngine<T> : ArrayRecordStore<T>.Writer
    {
        public WriterEngine(IGroup parentGroup):base(parentGroup)
        {
            parentGroup.InsertAttribute("DateCreated", DateTime.Now.ToString("o"));
            // FIXME
            parentGroup.InsertAttribute("ChangeListID", BuildVersion.VersionString);
        }
    }

    /// <summary>
    /// A base class to support the "PulseData" group format.  The format is characterized by
    /// field-parallel ragged arrays, concatenated across items (ZMWs) and indexed for item-specific access.
    /// Subordinate to the parent group containing the array fields, a "Markup" group holds scalar and fixed-width
    /// metadata corresponding to items.
    /// </summary>
    /// <typeparam name="T"></typeparam>
    /// <typeparam name="TRdr">Corresponding Reader class, containing method HdfIcdEntry attributes</typeparam>
    public class WriterBase<T,TRdr> : IDisposable where T:IZmwEvents
    {
        private static bool _compressData = true;

        /// <summary>
        /// Default setting whether to compress pulse and base datasets.
        /// For now we only compress on Linux -- for some reason compression doesn't work on windows
        /// </summary>
        public static bool CompressData
        {
            get { return _compressData; }
            set { _compressData = value; }
        }

        /// <summary>
        /// How much zlib compression to apply.  0 = none (fastest), 9 = max (slow)
        /// </summary>
        public static uint ZLibCompressionLevel = 3;

        // The parentGroup holds all parallel arrays, concatenated over items
        // ReSharper disable InconsistentNaming
        protected readonly IGroup parentGroup;

        // The markupGroup holds fixed-width, item-specific data
        protected readonly IGroup markupGroup;

        // The MetricGroup is optional, created on-demand
        private readonly string metricGroupName;
        private IGroup metricGroup;
        
        // The ArrayRecordStore writer
        protected readonly ArrayRecordStore<T>.Writer writer;
        // ReSharper restore InconsistentNaming

        // Private data used in constructing the writer
        private readonly List<string[]> markupContent;
        private readonly List<string[]> arraysContent;
        private readonly ArrayRecordStore<T>.Writer.ArrayGroup arraysWriter;

        /// <summary>
        /// Construct a Writer with the given group configuration.  Subclasses should use the Constructor Utils
        /// to define content/format.  The subclass constructor must call FinalizeContent() prior to completion.
        /// Default is to disable compression.
        /// </summary>
        /// <param name="fileGroup">The file group container</param>
        /// <param name="parentPath">The path of the Parent group, from the file root</param>
        /// <param name="markupGrpName">The name of the Markup group</param>
        /// <param name="metricGrpName">The name of the (optional) Metrics group</param>
        /// <param name="indexName">The name of the index dataset, stored under the Markup group</param>
        public WriterBase(IGroup fileGroup, string parentPath, string markupGrpName, string metricGrpName, string indexName) :
            this(fileGroup, parentPath, markupGrpName, metricGrpName, indexName, CompressData)
        {

        }

        /// <summary>
        /// Construct a Writer with the given group configuration.  Subclasses should use the Constructor Utils
        /// to define content/format.  The subclass constructor must call FinalizeContent() prior to completion.
        /// Default is to disable compression.
        /// </summary>
        /// <param name="fileGroup">The file group container</param>
        /// <param name="parentPath">The path of the Parent group, from the file root</param>
        /// <param name="markupGrpName">The name of the Markup group</param>
        /// <param name="metricGrpName">The name of the (optional) Metrics group</param>
        /// <param name="indexName">The name of the index dataset, stored under the Markup group</param>
        /// <param name="compression">Turn on compression gzip in the HDF5 file</param>
        public WriterBase(IGroup fileGroup, string parentPath, string markupGrpName, string metricGrpName,
                          string indexName, bool compression)
        {
            parentGroup = CreateParentGroup(fileGroup, parentPath);

            parentGroup.InsertAttribute("DateCreated", DateTime.Now.ToString("o"));
            // FIXME
            parentGroup.InsertAttribute("ChangeListID", BuildVersion.VersionString);

            // Get the file's dataset and attribute documentation
            readerDoc = InitReaderDoc();

            // Create the event-stream data writer
            writer = new ArrayRecordStore<T>.Writer(parentGroup, compression ? ZLibCompressionLevel : 0u);

            markupGroup = writer.CreateChildGroup(markupGrpName);
            metricGroupName = metricGrpName;

            markupContent = new List<string[]>();
            arraysContent = new List<string[]>();

            arraysWriter = writer.MakeParallelArrayGroup(indexName);
            markupContent.Add(new[] {indexName, MatType(typeof (int))});

            // Add standard Markup fields for ZMW identification
            AddMarkupSingleton<uint>("HoleNumber", b => (uint)b.Zmw.HoleNumber);
            AddMarkupFixedWidth<short>("HoleXY", b => new[] { (short)b.Zmw.X, (short)b.Zmw.Y });
            AddMarkupSingleton<Int16>("HoleChipLook", b => b.Zmw.HoleChipLook);

            var holeStatusDataset = AddMarkupSingleton<byte>("HoleStatus", b => (byte)b.Zmw.ZmwType);
            holeStatusDataset.WriteStringArrayAttribute("LookupTable", ZmwType.Types.Map(t => t.ToString()));
        }

        #region Utils for Writer Construcors

        private readonly Dictionary<string, HdfIcdEntry> readerDoc;

        // Initialize the ReaderDoc based on the Reader member attributes.
        private static Dictionary<string, HdfIcdEntry> InitReaderDoc()
        {
            try
            {
                // Readers should contain the static "GetIcd" method
                return (Dictionary<string, HdfIcdEntry>)(typeof(TRdr)).InvokeMember("GetIcd",
                                                                                     BindingFlags.InvokeMethod |
                                                                                     BindingFlags.Public |
                                                                                     BindingFlags.Static,
                                                                                     null, null, null);
            }
            catch
            {
                // If not, look for methods with HdfIcdEntry annotation
                var icd = new Dictionary<string, HdfIcdEntry>();

                var typeHdfIcdEntry = typeof(HdfIcdEntry);
                var readerMembers = (typeof(TRdr)).GetMembers().Where(
                    mi => mi.GetCustomAttributes(typeHdfIcdEntry, true).Length == 1).ToArray();

                foreach (var mi in readerMembers)
                {
                    var hdfEntry = (HdfIcdEntry)mi.GetCustomAttributes(typeHdfIcdEntry, true)[0];
                    var name = Path.GetFileName(hdfEntry.Path ?? mi.Name);

                    if (name == null)
                        continue;

                    icd.Add(name, hdfEntry);
                }

                return icd;
            }
        }

        protected HdfIcdEntry GetIcdEntry(string name, string description = null)
        {
            HdfIcdEntry icdEntry;
            if (readerDoc.TryGetValue(name, out icdEntry))
            {
                return icdEntry;
            }

            // For now, when no ICD entry is found, still allow the writer to provide description.
            if (description != null)
            {
                var icde = new HdfIcdEntry { Path = name, Detail = description };
                return icde;
            }

            throw new ApplicationException(String.Format("Field '{0}' is undocumented.", name));
        }

        protected IDataset AddMarkupSingleton<TR>(string name, Func<T, TR> getter)
        {
            // Force existence of documentation
            var icde = GetIcdEntry(name);

            var dataset = writer.WriteSetupSingleton(name, getter, markupGroup, icde.Detail, icde.Units);
            markupContent.Add(new[] { name, MatType(typeof(TR)) });
            return dataset;
        }

        protected void AddMarkupFixedWidth<TR>(string name, Func<T, TR[]> getter)
        {
            // Force existence of documentation
            var icde = GetIcdEntry(name);

            writer.WriteSetupFixedWidth(name, getter, markupGroup, icde.Detail, icde.Units);
            markupContent.Add(new[] { name, MatType(typeof(TR)) });
        }

        protected void AddArrayField<TR>(string name, Func<T, IList<TR>> getter, string description = null)
        {
            // Force existence of documentation
            var icde = GetIcdEntry(name, description);

            arraysWriter.AddArrayField(name, getter, null, icde.Detail, icde.Units);
            arraysContent.Add(new[] { name, MatType(typeof(TR)) });
        }

        protected void AddArrayField<TR>(string name, Func<T, TR[,]> getter, string description = null)
        {
            // Force existence of documentation
            var icde = GetIcdEntry(name, description);

            arraysWriter.AddArrayField(name, getter, null, icde.Detail, icde.Units);
            arraysContent.Add(new[] { name, MatType(typeof(TR)) });
        }

        protected void FinalizeContent()
        {
            // Finalize the ArraysWriter
            arraysWriter.Close(markupGroup);

            // Add Content description to the markupGroup
            string[,] content = MakeContentArray(markupContent);
            markupGroup.WriteStringArrayAttribute("Content", content);

            // Add Content description to the parentGroup
            content = MakeContentArray(arraysContent);
            parentGroup.WriteStringArrayAttribute("Content", content);
        }

        #endregion

        #region Public Methods

        /// <summary>
        /// Add an optional singleton ZMW metric
        /// </summary>
        public void AddMetricSingleton<TM>(string name, Func<T, TM> getter, string description = null)
        {
            InitMetricGroup();

            var icde = GetIcdEntry(name, description);
            writer.WriteSetupSingleton(name, getter, metricGroup, icde.Detail, icde.Units);
        }

        /// <summary>
        /// Add an optional fixed-width ZMW metric
        /// </summary>
        public void AddMetricFixedWidth<TM>(string name, Func<T, TM[]> getter, string description = null)
        {
            InitMetricGroup();

            var icde = GetIcdEntry(name, description);
            writer.WriteSetupFixedWidth(name, getter, metricGroup, icde.Detail, icde.Units);
        }

        public void AddMetricAnnotation<TA>(string key, TA value)
        {
            metricGroup.InsertAttribute(key, value);
        }

        /// <summary>
        /// Write a set of records
        /// </summary>
        /// <param name="data">The dataset to write to the file</param>
        /// <param name="chunkSize">The number of dataset elements to batch into a single HDF5 write.
        /// Set to zero to write the entire dataset in a single set of HDF5 writes
        /// </param>
        public void Write(IEnumerable<T> data, int chunkSize)
        {
            writer.WriteRecords(data,chunkSize);    
        }

        List<T> writeBuffer;

        /// <summary>
        /// For push-mode writing of data.  Give data one item at a time, it will be buffered and written.
        /// </summary>
        /// <param name="datum"></param>
        /// <param name="bufferSize"></param>
        public void BufferedWrite(T datum, int bufferSize)
        {
            if(writeBuffer == null)
                writeBuffer = new List<T>();

            writeBuffer.Add(datum);

            if (writeBuffer.Count >= bufferSize)
            {
                Write(writeBuffer, writeBuffer.Count);
                writeBuffer = null;
            }
        }

        private bool disposed = false;

        protected virtual void Dispose(bool disposing)
        {
            // Flush the write buffer if there's anything in it.
            if (writeBuffer != null && writeBuffer.Count > 0)
            {
                Write(writeBuffer, writeBuffer.Count);
                writeBuffer = null;
            }

            if (disposing && !disposed)
            {
                writer.Dispose();
                disposed = true;
            }
        }

        public void Dispose()
        {
            Dispose(true);
            GC.SuppressFinalize(this);
        }

        ~WriterBase()
        {
            Dispose(false);
        }

        #endregion

        #region Private Utils

        /// <summary>
        /// Create and return the parentGroup under the file group.
        /// The last group node is removed if it exists.
        /// </summary>
        /// <param name="fileGroup"></param>
        /// <param name="parentPath"></param>
        /// <returns></returns>
        private static IGroup CreateParentGroup(IGroup fileGroup, string parentPath)
        {
            var pgroup = fileGroup;
            var nodes = parentPath.Split('/');
            var iLast = nodes.Length - 1;

            for (int i = 0; i < iLast; i++)
                pgroup = pgroup.CreateGroup(nodes[i]);

            pgroup.RemoveChild(nodes[iLast]);
            return pgroup.CreateGroup(nodes[iLast]);
        }

        /// <summary>
        /// Util to re-format the content array
        /// </summary>
        /// <param name="contentList"></param>
        /// <returns></returns>
        private static string[,] MakeContentArray(IList<string[]> contentList)
        {
            int nItems = contentList.Count();
            var content = new string[2,nItems];

            for (int i=0; i<nItems; i++)
            {
                var v = contentList[i];
                content[0, i] = v[0];
                content[1, i] = v[1];
            }

            return content;
        }

        /// <summary>
        /// Initialize the optional MetricGroup
        /// </summary>
        private void InitMetricGroup()
        {
            if (metricGroup == null)
                metricGroup = writer.CreateChildGroup(metricGroupName);
        }

        private static string MatType(Type type)
        {
            var tstr = type.ToString();
            var t = tstr.StartsWith("System.") ? tstr.Substring(7).ToLower() : tstr;

            if (t.Equals("byte"))
                t = "uint8";

            return t;
        }

        #endregion
    }

    /// <summary>
    /// Writes a stream of IZmwConsensusBases objects to HDF5, along with required metadata.
    /// </summary>
    public class ConsensusBaseWriter : WriterBase<IZmwConsensusBases, ConsensusBaseReader>
    {
        // Private data used in constructing the writer
        private readonly IGroup passesGroup;
        private readonly ArrayRecordStore<IZmwConsensusBases>.Writer.ArrayGroup passesWriter;

        public ConsensusBaseWriter(IGroup fileGroup, string groupName)
            : base(fileGroup, "PulseData/" + groupName, "ZMW", "ZMWMetrics", "NumEvent")
        {
            // Group attributes
            parentGroup.InsertAttribute("SchemaRevision", "1.0");
            parentGroup.InsertAttribute("QVDecoding", "QV = -10 * log10(p), where p is the probability of error");

            // Array fields
            AddArrayField<byte>("Basecall", b => b.Base.Map(v => (byte) v));
            AddArrayField<byte>("QualityValue", b => b.QV);
            AddArrayField<byte>("InsertionQV", b => b.InsertionQV);
            AddArrayField<byte>("DeletionQV", b => b.DeletionQV);
            AddArrayField<byte>("DeletionTag", b => b.DeletionTag.Map(v => (byte) v));
            AddArrayField<byte>("SubstitutionQV", b => b.SubstitutionQV);
            AddArrayField<byte>("SubstitutionTag", b => b.SubstitutionTag.Map(v => (byte) v));

            // Passes group
            passesGroup = writer.CreateChildGroup("Passes");
            passesWriter = writer.MakeParallelArrayGroup("NumPasses");

            AddPassesField<uint>("PassStartBase", b => b.InsertRegions.Map(r => (uint)r.Start));
            AddPassesField<uint>("PassNumBases", b => b.InsertRegions.Map(r => (uint)r.Length));
            AddPassesField<byte>("PassDirection", b => b.InsertRegions.Map(r => (byte)r.Strand));
            AddPassesField<byte>("AdapterHitBefore", b => b.InsertRegions.Map(r => r.AdapterHitBefore ? (byte) 1 : (byte) 0));
            AddPassesField<byte>("AdapterHitAfter",
                                 b => b.InsertRegions.Map(r => r.AdapterHitAfter ? (byte) 1 : (byte) 0));

            passesWriter.Close(passesGroup);
            
            // Write the predicted accuracy and insert length
            //AddMetricSingleton<float>("PredictedAccuracy", b => b.PredictedAccuracy);
            //AddMetricSingleton<uint>("InsertReadLength", b => (uint) b.InsertReadLength);
            AddMetricSingleton<byte>("Productivity", b => (byte)b.Metrics.Productivity);
            AddMetricSingleton<float>("ReadScore", b => b.Metrics.ReadScore);


            // Write the regions table to the bases -- goes in /PulseData/Regions
            if (parentGroup != null)
            {
                var ww = new RegionAnnotator.Writer(fileGroup.CreateGroup("PulseData"));
                Action<IZmwBases> writeRegions = bas => ww.Write(bas.Metrics.Regions);
                writer.AddExtraWriter(bases => bases.ForEach(writeRegions));
            }

            // Finish building the writer
            FinalizeContent();
        }

        private void DeleteAttribute(string path, string name)
        {
            var grp = (IGroup)parentGroup.File;

            foreach (var node in path.Split(new[]{'/'}))
            {
                grp = (IGroup)grp.GetChild(node);

                if (grp == null)
                    return;
            }

            grp.DeleteAttribute(name);
        }

        private void WriteAttribute(string path, string name, object obj)
        {
            var grp = (IGroup)parentGroup.File;

            foreach (var node in path.Split(new[]{'/'}))
            {
                grp = grp.CreateGroup(node);
            }

            grp.InsertAttribute(name, obj);
        }

        public void AddChemistryInformation(string sequencingChemistry)
        {
            DeleteChemistryInformation();
            WriteAttribute("ScanData/RunInfo", "SequencingChemistry", sequencingChemistry);
        }

        public void AddChemistryInformation(string bindingKit, string sequencingKit, string changeListId)
        {
            DeleteChemistryInformation();
            WriteAttribute("ScanData/RunInfo", "BindingKit", bindingKit);
            WriteAttribute("ScanData/RunInfo", "SequencingKit", sequencingKit);
            WriteAttribute("PulseData/BaseCalls", "ChangeListID", changeListId);
        }

        private void DeleteChemistryInformation()
        {
            DeleteAttribute("ScanData/RunInfo", "SequencingChemistry");
            DeleteAttribute("ScanData/RunInfo", "BindingKit");
            DeleteAttribute("ScanData/RunInfo", "SequencingKit");
            DeleteAttribute("PulseData/BaseCalls", "ChangeListID");
        }

        private void AddPassesField<TR>(string name, Func<IZmwConsensusBases, IList<TR>> getter)
        {
            // Force existence of documentation
            var icde = GetIcdEntry(name);
            passesWriter.AddArrayField(name, getter, passesGroup, icde.Detail, icde.Units);
        }
    }

    /// <summary>
    /// Writes a stream of IZmwBase objects to HDF5, along with required metadata.
    /// </summary>
    public class BaseWriter : WriterBase<IZmwBases, BaseReader>
    {
        /// <summary>
        /// A static method to construct and configure a BaseWriter for standard PulseToBase output.
        /// For behavior consistent with latest P2B, implementations should use this method rather than
        /// "new BaseWriter(...)".
        /// </summary>
        /// <param name="baseFile"></param>
        /// <param name="regionFile"></param>
        /// <param name="metadata"></param>
        /// <param name="richQvs"></param>
        /// <param name="kineticData"></param>
        /// <returns></returns>
        public static BaseWriter ForPulseToBase(IChunkFile baseFile, IChunkFile regionFile, IMovieMetadata metadata,
                                                bool richQvs = true, bool kineticData = true)
        {
            // Construct the base writer
            var bw = new BaseWriter(baseFile, regionFile, richQvs, kineticData);

            // Configure with optional ZMW metrics
            AddBasecallsMetrics(bw, metadata, richQvs);

            return bw;
        }

        public static bool WriteTdMetrics = true;

        public static void AddBasecallsMetrics(BaseWriter bw, IMovieMetadata md, bool writeRichQVs)
        {
            // Conversion for signal-level metrics into counts.
            float conv = 1.0f;
            Func<float[], float[]> toCounts = p => p.Map(v => (conv*v));

            bw.AddMetricFixedWidth<float>("BaseFraction", b => b.Metrics.BaseFraction);
            bw.AddMetricSingleton<float>("BaseRate", b => b.Metrics.BaseRate);
            bw.AddMetricSingleton<float>("BaseWidth", b => b.Metrics.BaseWidth);
            bw.AddMetricSingleton<float>("BaseIpd", b => b.Metrics.BaseIpd);
            bw.AddMetricSingleton<float>("Pausiness", b => b.Metrics.Pausiness);
            bw.AddMetricSingleton<float>("LocalBaseRate", b => b.Metrics.LocalBaseRate);
            bw.AddMetricSingleton<float>("DarkBaseRate", b => b.Metrics.DarkBaseRate);
            bw.AddMetricSingleton<byte>("Productivity", b => (byte) b.Metrics.Productivity);
            bw.AddMetricSingleton<byte>("Loading", b => (byte) b.Metrics.Loading);
            bw.AddMetricSingleton<float>("ReadScore", b => b.Metrics.ReadScore);
            bw.AddMetricFixedWidth<float>("HQRegionSNR", b => b.Metrics.HQRegionSNR);
            bw.AddMetricFixedWidth<float>("HQRegionIntraPulseStd", b => toCounts(b.Metrics.HQRegionIntraPulseStd));
            bw.AddMetricSingleton<float>("HQRegionStartTime", b => b.Metrics.HQRegionStartTime);
            bw.AddMetricSingleton<float>("HQRegionEndTime", b => b.Metrics.HQRegionEndTime);
            bw.AddMetricFixedWidth<float>("HQRegionEstPkmid", b => toCounts(b.Metrics.HQRegionEstPkmid));
            bw.AddMetricFixedWidth<float>("HQRegionEstPkstd", b => toCounts(b.Metrics.HQRegionEstPkstd));
            bw.AddMetricFixedWidth<float>("HQRegionPkzvar", b => b.Metrics.HQRegionPkzvar);
            bw.AddMetricSingleton<float>("RmBasQv", b => b.Metrics.RmBasQv);
            bw.AddMetricFixedWidth<float>("CmBasQv", b => b.Metrics.CmBasQv);

            if (writeRichQVs)
            {
                bw.AddMetricFixedWidth<float>("CmInsQv", b => b.Metrics.CmInsQv);
                bw.AddMetricFixedWidth<float>("CmDelQv", b => b.Metrics.CmDelQv);
                bw.AddMetricFixedWidth<float>("CmSubQv", b => b.Metrics.CmSubQv);

                bw.AddMetricSingleton<float>("RmInsQv", b => b.Metrics.RmInsQv);
                bw.AddMetricSingleton<float>("RmDelQv", b => b.Metrics.RmDelQv);
                bw.AddMetricSingleton<float>("RmSubQv", b => b.Metrics.RmSubQv);
            }

            if (WriteTdMetrics)
            {
                bw.AddMetricFixedWidth<short>("NumBaseVsT", b => b.Metrics.NumBaseVsT);
                bw.AddMetricFixedWidth<short>("NumPauseVsT", b => b.Metrics.NumPauseVsT);
                bw.AddMetricFixedWidth<float>("BaseRateVsT", b => b.Metrics.BaseRateVsT);

                bw.AddMetricAnnotation("TdmIntervalSec", ZmwMetricsBases.CfgTdmIntervalSec);
                bw.AddMetricAnnotation("PauseIpdThreshSec", ZmwMetricsBases.CfgPauseIpdThreshSec);
            }
        }

        public BaseWriter(IGroup fileGroup, bool writeRichQvs, bool writeKineticData) :
            this(fileGroup, null, writeRichQvs, writeKineticData)
        {

        }

        /// <summary>
        /// Setup a Basecalls writer. A group called PulseData/BaseCalls will be created below the HDF5 group passed
        /// as to the constructor.  The basecalls datasets will be stored in that group, according to the 
        /// Base Calls ICD. When using this writer it is important to close the file to ensure the data is flushed to disk.
        /// </summary>
        /// <param name="fileGroup">The group beneath which to create the PulseData/BaseCalls output group. For a standard 
        ///   basecalls HDF5 file this will be the top-level file group.</param>
        /// <param name="regionFile"></param>
        /// <param name="writeRichQvs">Whether or not to write rich base QVs</param>
        /// <param name="writeKineticData">Whether or not to write kinetic data</param>
        public BaseWriter(IGroup fileGroup, IChunkFile regionFile, bool writeRichQvs, bool writeKineticData)
            : base(fileGroup, "PulseData/BaseCalls", "ZMW", "ZMWMetrics", "NumEvent")
        {
            // Group attributes            
            parentGroup.InsertAttribute("SchemaRevision", "1.0");
            parentGroup.InsertAttribute("QVDecoding",
                                        "Standard Phred encoding: QV = -10 * log10(p) - where p is the probability of error");

            // Array fields
            AddArrayField<byte>("Basecall", b => b.Base.Map(v => (byte) v));
            AddArrayField<int>("PulseIndex", b => b.PulseIndex);
            AddArrayField<byte>("QualityValue", b => b.QV, "Called base quality value");

            if (writeRichQvs)
            {
                AddArrayField<byte>("InsertionQV", b => b.InsertionQV);
                AddArrayField<byte>("DeletionQV", b => b.DeletionQV);
                AddArrayField<byte>("DeletionTag", b => b.DeletionTag.Map(v => (byte) v));
                AddArrayField<byte>("SubstitutionQV", b => b.SubstitutionQV);
                AddArrayField<byte>("SubstitutionTag", b => b.SubstitutionTag.Map(v => (byte) v));
                AddArrayField<byte>("MergeQV", b => b.MergeQV);
            }

            if (writeKineticData)
            {
                AddArrayField<ushort>("WidthInFrames", b => b.WidthInFrames);
                AddArrayField<ushort>("PreBaseFrames", b => b.PreBaseFrames);
            }

            // Write the regions table to the bases, if a regionGroup was supplied
            if (regionFile != null)
            {
                var ww = new RegionAnnotator.Writer(regionFile);
                Action<IZmwBases> writeRegions = bas => ww.Write(bas.Metrics.Regions);
                writer.AddExtraWriter(bases => bases.ForEach(writeRegions));
            }

            // Finish building the writer
            FinalizeContent();
        }
    }
}
