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
using System.IO;
using System.Text.RegularExpressions;
using PacBio.IO.ICD;
using PacBio.HDF;
using PacBio.Utils;

namespace PacBio.IO
{
    /// <summary>
    /// A base class for constructing readers for per-ZMW event stream HDF5 files,
    /// (e.g. pls.h5, bas.h5, etc) where the source consists of a single file. 
    /// </summary>
    /// <typeparam name="TOutputType">The type of object being read</typeparam>
    public abstract class ReaderBase<TOutputType> : DataSource<TOutputType> where TOutputType : class
    {
        #region Members

        // Some metadata about the contents of the file
        private string changelist;
        private DateTime dateCreated;
        protected float SchemaRevision;

        // Given an index into the the ZMW array, generate a thunk that returns the output object
        protected Func<int, Func<TOutputType>> ZmwEventFunc;

        #endregion

        #region Structors
        // ReSharper disable DoNotCallOverridableMethodsInConstructor

        protected ReaderBase(string filename, string dataGroupName)
        {
            FileName = filename;

            // Use a reasonable chunk cache of 5MB when opening, the defaults may be too small.
            var fileAccessProps = new HDFFileAccessProperty();
            fileAccessProps.SetCache(10000, 1023, (ulong) 5e6, 0.75);

            // Open file read-only
            Chunks = HDFFile.Open(filename, FileMode.Open, FileAccess.Read, fileAccessProps);
            
            Init(dataGroupName);
        }

        protected ReaderBase(IChunkFile chunks, string filename, string dataGroupName)
        {
            FileName = filename;
            Chunks = chunks;
            Init(dataGroupName);
        }

        private void Init(string dataGroupName)
        {
            try
            {
                // Get the HDF groups that we want
                try
                {
                    var pulseData = (IGroup) Chunks.GetChild("PulseData");
                    EventGroup = (IGroup)pulseData.GetChild(dataGroupName);
                    ZmwGroup = (IGroup)EventGroup.GetChild("ZMW");

                    // Optional ZMW Metrics group
                    ZmwOptGroup = EventGroup.GetChild("ZMWMetrics") as IGroup;
                }
                catch (NullReferenceException)
                {
                    var msg = String.Format("HDF data file {0} appears to be corrupted.", FileName);
                    throw new InvalidDataException(msg);
                }

                // Use the passed in ZmwSource, or load it from the current file if it's null.
                if (Chunks.GetChild("ScanData") != null)
                    ZmwSource = MetadataReader.GetMetadataReader((IGroup) Chunks.GetChild("ScanData"), ZmwGroup);
                else
                    throw new ArgumentException(
                        "No IZmwSource was passed, and ScanData group wasn't available in current file.");

                // Software version / changeList
                var clAttr = EventGroup.GetAttribute("ChangeListID");
                changelist = clAttr != null ? clAttr.ReadSingleton<string>() : "Unknown version";

                dateCreated = Helpers.GetDateAttribute(EventGroup, "DateCreated");

                // The 'SchemaRevision' attribute was introduced in 1.2.2 it is the official indicator of the
                // layout of the dataset.  Prior to 1.2.2, the 'Version' attribute did the same thing.
                var schemaRevisionAttr = EventGroup.GetAttribute("SchemaRevision") ?? EventGroup.GetAttribute("Version");
                SchemaRevision = float.Parse(schemaRevisionAttr != null ? schemaRevisionAttr.ReadSingleton<string>() : "0");

                // Pull out the NumPulses data
                var numPulseDataset = (IDataset) ZmwGroup.GetChild("NumEvent");
                var a = numPulseDataset.Read();

                // ReSharper disable PossibleNullReferenceException
                ZmwNumEvents = a is Int32[]
                    ? ((Int32[])a).Select(i => (uint) i).ToArray()
                    : (UInt32[]) a;
                // ReSharper restore PossibleNullReferenceException

                NumZmws = ZmwNumEvents.Length;

                // Fill out an array containing the start offset of all the zmws
                ZmwStartIndex = new long[NumZmws];
                long startIdx = 0;
                for (int i = 0; i < ZmwNumEvents.Length; i++)
                {
                    ZmwStartIndex[i] = startIdx;
                    startIdx += ZmwNumEvents[i];
                }
            }
            catch (IOException e)
            {
                throw new IOException(String.Format("Error Reading HDF Pulse file: {0}", FileName), e);
            }
        }

        protected override void Dispose(bool disposing)
        {
            if (disposing && Chunks != null)
            {
                Chunks.Dispose();
                Chunks = null;
            }
        }

        // ReSharper restore DoNotCallOverridableMethodsInConstructor
        #endregion

        #region Properties

        internal IChunkFile Chunks { get; private set; }
        private IGroup ZmwGroup { get; set; }
        private IGroup ZmwOptGroup { get; set; }
        protected IGroup EventGroup { get; set; }
        protected int NumZmws { get; set; }
        private long[] ZmwStartIndex { get; set; }

        public double FrameRate { get { return ZmwSource.FrameRate; } }

        public string[] Fields
        {
            get
            {
                return EventGroup.GetChildren().Select(ce => ce.Name.Split('/').Last()).ToArray(); 
            }
        }

        // Direct access to the ZMW Markup data is provided here, even though in practice
        // only the ZmwNumEvents property has been used.  The others are convenience methods,
        // and provide the annotation point for the ICD documentation.  The NumEvent field
        // is special, because it is the IndexField for data in the Event group.  The Description
        // for that data cannot be provided by this Reader because it is written differently than
        // the others.

        [HdfIcdEntry(Path = "ZMW/HoleNumber", Detail = "Number assigned to each ZMW on the chip")]
        public int[] HoleNumber { get { return ZmwSource.ZmwIndexer.HoleNumber; } }

        [HdfIcdEntry(Path = "ZMW/HoleXY", Detail = "Grid coordinates assigned to each ZMW on the chip")]
        public XYPoint[] HoleXY { get { return ZmwSource.ZmwIndexer.HoleXY; } }

        [HdfIcdEntry(Path = "ZMW/HoleChipLook", Detail = "Look of ZMW in chip layout")]
        public Int16[] HoleChipLook { get { return ZmwSource.ZmwIndexer.HoleChipLook; } }

        [HdfIcdEntry(Path = "ZMW/HoleStatus", Detail = "Type of ZMW that produced the data")]
        public ZmwType[] HoleStatus { get { return ZmwSource.ZmwIndexer.ZmwType; } }

        [HdfIcdEntry(Path = "ZMW/NumEvent", Detail = "Event counts per ZMW for all fields in the Event group")]
        public uint[] ZmwNumEvents { get; set; }
        
        #endregion

        #region Utility Methods

        /// <summary>
        /// Generate a custom 'serialization' of an Enum type.
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <returns></returns>
        public static string[] MappingForEnum<T>()
        {
            var names = Enum.GetNames(typeof(T));
            return names.Select(n => String.Format("{0}:{1}", (int)Enum.Parse(typeof(T), n), n)).ToArray();
        }

        /// <summary>
        /// Given a mapping array, with each element of the form "0:Label",
        /// return a function that maps byte values to enum values, as long as
        /// "Label" is a defined value for the enum type.
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="mapping"></param>
        /// <param name="undefinedVal"></param>
        /// <returns></returns>
        public static Func<byte, T> ByteToEnum<T>(string[] mapping, T undefinedVal)
        {
            // Build the dictionary
            var reng = new Regex(@"([0-9]+):([a-z_][a-z0-9_]*)", RegexOptions.IgnoreCase);
            var dict = new Dictionary<byte, T>();

            foreach (var expr in mapping)
            {
                var matches = reng.Match(expr);
                if (!matches.Success)
                    throw new ApplicationException("Invalid Enum mapping: bad expression.");

                var key = (byte)Int32.Parse(matches.Groups[1].Value);
                T value;
                try
                {
                    value = (T)Enum.Parse(typeof(T), matches.Groups[2].Value);
                }
                catch (ArgumentException)
                {
                    value = undefinedVal;
                }

                if (dict.ContainsKey(key))
                    throw new ApplicationException("Invalid Enum mapping: duplicate entry.");

                dict.Add(key, value);
            }

            // Return a closure
            return v => dict.ContainsKey(v) ? dict[v] : undefinedVal;
        }

        public uint GetNumEvents(int idx)
        {
            return ZmwNumEvents[idx];
        }

        private void CheckType(Type targetType, Type nativeType)
        {
            if (targetType != nativeType)
            {
                string msg = String.Format("Requested return type '{0}' does not match stored type '{1}'",
                                           targetType.FullName, nativeType.FullName);
                throw new ArgumentException(msg);
            }
        }
        
        /// <summary>
        /// Build a function that will return the feature data for a zmw event.  The returned function takes 
        /// the holeNumber as an argument, and returns the data for the single zmw. Use this method to access
        /// pulse features with a scalar value per pulse
        /// </summary>
        /// <typeparam name="T">The return type - must match the type stored in the HDF5 file</typeparam>
        /// <param name="pulseFeature">The name of the event feature in the HDF5 file</param>
        protected Func<int, T[]> MakeFeatureLookup1D<T>(string pulseFeature)
        {
            var dataset = (IDataset)EventGroup.GetChild(pulseFeature);

            if (dataset == null)
            {
                string msg = String.Format("HDF file does not contain the requested field: {0}", pulseFeature);
                throw new ArgumentException(msg);
            }

            var nativeType = dataset.Datatype.NativeType;

            CheckType(typeof(T), nativeType);
            
            var ones2D = new long[] { 1, 1 };
            var ones1D = new long[] { 1 };

            var dataspace = dataset.Dataspace;
            if ((dataspace.Dimensions.Length > 1 && dataspace.Dimensions[1] != 1) || dataspace.Dimensions.Length > 2)
            {
                throw new ArgumentException(String.Format("HDF array '{0}' is not 1D as required", pulseFeature));
            }

            Func<int, T[]> res;

            if (dataspace.Dimensions.Length == 2)
            {
                // A closure that pulls out the data for the given ZMW index, for the given pulseFeature
                res = delegate(int index0)
                          {
                              var start = new [] { ZmwStartIndex[index0], 0 };
                              var count = new long[] { ZmwNumEvents[index0], 1 };
                              var outCount = new long[] { ZmwNumEvents[index0] };
                              using (var rspace = dataset.Dataspace)
                              {
                                  // Make a new dataspace to be used for this access
                                  rspace.SelectHyperslab(start, ones2D, count, ones2D);

                                  var target = Array.CreateInstance(nativeType, outCount);
                                  dataset.Read(ref target, rspace);
                                  return (T[]) target;
                              }
                      };
            }
            else
            {
                // A closure that pulls out the data for the given ZMW index, for the given pulseFeature
                res = delegate(int index0)
                          {
                              var start = new [] { ZmwStartIndex[index0] };
                              var count = new long[] { ZmwNumEvents[index0] };
                              using (var rspace = dataset.Dataspace)
                              {
                                  rspace.SelectHyperslab(start, ones1D, count, ones1D);

                                  var target = Array.CreateInstance(nativeType, count);
                                  dataset.Read(ref target, rspace);
                                  return (T[]) target;
                              }
                      };
            }

            return res;
        }

        /// <summary>
        /// Build a function that will return the feature data for a zmw event.  The returned function takes 
        /// the holeNumber as an argument, and returns the data for the single zmw. Use this method to access
        /// pulse features with a scalar value per pulse.  Use this method when files with a legacy type exist;
        /// when those are encountered, the data will be read and cast to the (new) target type.
        /// </summary>
        /// <typeparam name="T">The return type</typeparam>
        /// <typeparam name="TLgcy">The type contained in the HDF5 file</typeparam>
        /// <param name="pulseFeature">The name of the event feature in the HDF5 file</param>
        /// <param name="convert">A funtion to convert values from TLgcy to T</param>
        protected Func<int, T[]> MakeFeatureLookup1D<T, TLgcy>(string pulseFeature, Func<TLgcy,T> convert)
        {
            var dataset = (IDataset)EventGroup.GetChild(pulseFeature);

            if (dataset == null)
            {
                string msg = String.Format("HDF file does not contain the requested field: {0}", pulseFeature);
                throw new ArgumentException(msg);
            }

            var nativeType = dataset.Datatype.NativeType;
            //var targetType = typeof (T);
            var legacyType = typeof (TLgcy);

            if (nativeType == legacyType)
            {
                // Return a closure that reads as legacy type and casts to target type
                Func<int, T[]> res = index0 => 
                          {
                              var f = MakeFeatureLookup1D<TLgcy>(pulseFeature);
                              var tLgcy = f(index0);
                              return tLgcy.Map(convert);
                          };

                return res;
            }

            return MakeFeatureLookup1D<T>(pulseFeature);
        }

        /// <summary>
        /// Build a function that will return the feature data for a zmw event.  The returned function takes 
        /// the holeNumber as an argument, and returns the data for the single zmw. Use this method to access
        /// pulse features, where the feature itself is a fixed-length vector
        /// </summary>
        /// <typeparam name="T">The HDF5 file datatype</typeparam>
        /// <typeparam name="TR">The return HDF5 dataset name</typeparam>
        /// <param name="pulseFeature"></param>
        /// <param name="f"></param>
        /// <returns></returns>
        protected Func<int, TR[,]> MakeFeatureLookup2D<T, TR>(string pulseFeature, Func<T, TR> f)
        {
            var dataset = (IDataset)EventGroup.GetChild(pulseFeature);

            if (dataset == null)
            {
                string msg = String.Format("HDF file does not contain the requested field: {0}", pulseFeature);
                throw new ArgumentException(msg);
            }

            var nativeType = dataset.Datatype.NativeType;

            CheckType(typeof(T), nativeType);

            var dataspace = dataset.Dataspace;
            if (dataspace.Dimensions.Length != 2)
            {
                throw new ArgumentException(String.Format("HDF array '{0}' is not 2D as required", pulseFeature));
            }

            var ones = new long[] {1, 1};
            var width = dataset.Dataspace.Dimensions[1];

            Func<int, TR[,]> res = delegate(int index0)
                                       {
                                           var start = new[] {ZmwStartIndex[index0], 0};
                                           var count = new[] {ZmwNumEvents[index0], width};
                                           using (var rspace = dataset.Dataspace)
                                           {
                                               rspace.SelectHyperslab(start, ones, count, ones);

                                               var target = Array.CreateInstance(nativeType, count);
                                               dataset.Read(ref target, rspace);

                                               var intermediate = target as T[,];
                                               return intermediate.Map(f);
                                           }
                                       };

            return res;
        }

        /// <summary>
        /// Build a function to lookup a vector zmw metric, given a holeNumber
        /// </summary>
        /// <typeparam name="T">The HDF5 dataset type</typeparam>
        /// <typeparam name="TR">The return vector type</typeparam>
        /// <param name="traceMetric">The name of the trace metric</param>
        /// <param name="f">A function to convert from the file type to the return type</param>
        /// <returns></returns>
        protected Func<int, TR[]> MakeMetricLookup2D<T,TR>(string traceMetric, Func<T, TR> f)
        {
            var dataset = (IDataset)ZmwGroup.GetChild(traceMetric);

            // Check the dataset type
            var nativeType = dataset.Datatype.NativeType;
            CheckType(typeof(T), nativeType);

            // Check the dataset dimensions
            var dataspace = dataset.Dataspace;
            if (dataspace.Dimensions.Length != 2)
            {
                throw new ArgumentException(String.Format("HDF array '{0}' is not 2D as required", traceMetric));
            }

            var length = dataset.Dataspace.Dimensions[0];
            var width = dataset.Dataspace.Dimensions[1];
            T[,] data = null;
            // A closure that obtains the data for the given ZMW index, using weak ref for the full dataset
            return delegate(int index0)
                       {                           
                           if (data == null)
                           {
                               data = new T[length,width];
                               Array outref = data;
                               dataset.Read(ref outref, dataspace);
                           }

                           return ((int) width).Fill(i => data[index0, i]).Map(f);
                       };
        }

        /// <summary>
        /// Get an optional ZMWMetric dataset
        /// </summary>
        /// <param name="metricName"></param>
        /// <returns></returns>
        private IDataset OptZmwMetricDataset(string metricName)
        {
            return ZmwOptGroup == null ? null : ZmwOptGroup.GetChild(metricName) as IDataset;
        }

        /// <summary>
        /// Read an optional metric annotation, or return the default value
        /// if the data is missing.
        /// </summary>
        /// <typeparam name="TR"></typeparam>
        /// <param name="key"></param>
        /// <param name="defValue"></param>
        /// <returns></returns>
        protected TR ReadOptMetricAnnotation<TR>(string key, TR defValue)
        {
            if (ZmwOptGroup == null)
                return defValue;

            var attr = ZmwOptGroup.GetAttribute(key);
            return (attr != null ? attr.ReadSingleton<TR>() : defValue);
        }

        // Initialize a lookup table for a classification metric
        protected Func<byte,TEnum> DefineClassificationMapping<TEnum>(string metricName, TEnum undefinedVal)                                                                     
        {
            var dataset = OptZmwMetricDataset(metricName);
            if (dataset == null)
                return null;
            
            // For this type of metric (a classification), this attribute should be an array
            // of strings, with each element in the format "0:Label"
            var attr = dataset.GetAttribute("UnitsOrEncoding");
            if (attr == null)
                return null;

            // In this case, this attribute should be a string array.
            var mapping = attr.Read() as string[];

            return mapping == null ? null : ByteToEnum(mapping, undefinedVal);
        }

       
        /// <summary>
        /// Build a function that will return an optional singleton metric for a ZMW.
        /// The returned function takes the hole number as an argument and returns the metric.
        /// </summary>
        /// <typeparam name="TR">The datatype of the returned metric (must the the same as stored in HDF5)</typeparam>
        /// <param name="metricName">The HDF5 dataset name for the metric, stored optionally under .../ZMWMetrics</param>
        /// <param name="f">A function to compute the metric on-the-fly when it does not exist in the HDF5 file</param>
        /// <returns></returns>
        protected Func<int, TR> MakeOptMetricLookupSingleton<TR>(string metricName, Func<TOutputType, TR> f)
        {
            var dataset = OptZmwMetricDataset(metricName);
            if (dataset == null)
            {
                // Return the function that will compute on-the-fly
                return index0 => f(ZmwEventFunc(index0)());
            }

            // Check the dataset type
            var nativeType = dataset.Datatype.NativeType;
            CheckType(typeof(TR), nativeType);

            // Check the dataset dimensions
            IDataspace dataspace = dataset.Dataspace;
            if (dataspace.Dimensions.Length != 1)
            {
                throw new ArgumentException(String.Format("HDF array '{0}' is not 1D as required", metricName));
            }

            long length = dataset.Dataspace.Dimensions[0];
            TR[] data = null;
            // A closure that obtains the data for the given ZMW index, using weak ref for the full dataset
            return delegate(int index0)
                       {                          
                           if (data == null)
                           {
                               data = new TR[length];
                               Array outref = data;
                               dataset.Read(ref outref, dataspace);
                           }

                           return data[index0];
                       };
        }

        /// <summary>
        /// Build a function that will return an optional fixed-width metric for a ZMW.
        /// The returned function takes the hole number as an argument and returns the metric.
        /// </summary>
        /// <typeparam name="TR">The datatype of the returned metric (must the the same as stored in HDF5)</typeparam>
        /// <param name="metricName">The HDF5 dataset name for the metric, stored optionally under .../ZMWMetrics</param>
        /// <param name="f">A function to compute the metric on-the-fly when it does not exist in the HDF5 file</param>
        /// <returns></returns>
        protected Func<int, TR[]> MakeOptMetricLookupFixedWidth<TR>(string metricName, Func<TOutputType, TR[]> f)
        {
            var dataset = OptZmwMetricDataset(metricName);
            if (dataset == null)
            {
                // Return a function that will compute on-the-fly
                return index0 => f(ZmwEventFunc(index0)());
            }

            // Check the dataset type
            var nativeType = dataset.Datatype.NativeType;
            CheckType(typeof (TR), nativeType);

            // Check the dataset dimensions
            IDataspace dataspace = dataset.Dataspace;
            if (dataspace.Dimensions.Length != 2)
            {
                throw new ArgumentException(String.Format(
                    "HDF array '{0}' is not 2D as required", metricName));
            }

            long length = dataset.Dataspace.Dimensions[0];
            long width = dataset.Dataspace.Dimensions[1];

            // A closure that obtains the data for the given ZMW index, using weak ref for the full dataset
            TR[,] data =null;
            return delegate(int index0)
                       {
                           
                           if (data == null)
                           {
                               data = new TR[length,width];
                               Array outref = data;
                               dataset.Read(ref outref, dataspace);
                           }

                           return ((int) width).Fill(i => data[index0, i]);
                       };
        }

        #endregion

        #region IList

        public override TOutputType this[int index]
        {
            get
            {
                //return zmwEventArray[index]();
                return ZmwEventFunc(index)();
            }
            set
            {
                throw new NotImplementedException();
            }
        }

        #endregion

        #region ISourceIdentifier

        [HdfIcdEntry(Path = "ChangeListID", Detail="Revision ID of software that produced this data",
                     Units="<Major>.<Minor>.<Micro>.<Hotfix>.<changeNumber>")]
        public override string ChangelistID
        {
            get { return changelist; }
        }

        public override string SoftwareVersion
        {
            get { return changelist; }
        }

        public override string FileName { get; protected set;}

        [HdfIcdEntry(Detail = "Creation date-time of the data group", Units = "ISO 8601")]
        public override DateTime DateCreated
        {
            get { return dateCreated; }
        }

        #endregion
    }
}
