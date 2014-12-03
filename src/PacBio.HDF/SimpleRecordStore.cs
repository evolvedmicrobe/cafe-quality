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
using System.Reflection;
using PacBio.Utils;

namespace PacBio.HDF
{
    public class SingletonReader<R>
    {
        public static long[] la(int i)
        {
            return new long[1] { i };
        }

        private IDataset dataset;

        public SingletonReader(IDataset ds)
        {
            dataset = ds;
        }

        public R[] Read(int start, int n)
        {
            using (var dspace = dataset.Dataspace)
            {
                dspace.SelectHyperslab(la(start), la(1), la(n), null);

                R[] indata = new R[n];
                Array outref = indata;
                dataset.Read(ref outref, dspace);

                return indata;
            }
        }
    }

    public class ArrayReader<R>
    {
        public static long[] la(int i)
        {
            return new long[1] { i };
        }

        private IDataset dataset;
        private int[] recordSizes;

        public ArrayReader(IDataset ds, int[] rs)
        {
            dataset = ds;
            recordSizes = rs;
        }

        public R[][] Read(int start, int n)
        {
            int startPoint = 0;
            for (int i = 0; i < start; i++)
                startPoint += recordSizes[i];

            int[] sizes = new int[n]; //recordSizes.Slice(start, n);
            Array.Copy(recordSizes, start, sizes, 0, n);

            int readSize = sizes.Sum();
            var startPoints = new int[sizes.Length];
            startPoints[0] = 0;
            for (int i = 1; i < startPoints.Length; i++)
            {
                startPoints[i] = startPoints[i - 1] + sizes[i - 1];
            }

            R[] readData = new R[readSize];
            Array readRef = readData;
            using (var dspace = dataset.Dataspace)
            {
                dspace.SelectHyperslab(la(startPoint), la(1), la(readSize), null);
                dataset.Read(ref readRef, dspace);

                var result = new R[startPoints.Length][];

                for (int i = 0; i < startPoints.Length; i++)
                {
                    var rr = new R[sizes[i]];
                    Array.Copy(readData, startPoints[i], rr, 0, sizes[i]);
                    result[i] = rr;
                }

                return result;
            }
        }
    }

    /// <summary>
    /// A class for storing simple C# types to parallel HDF5 arrays.
    /// </summary>
    public class SimpleRecordStore
    {
        public static readonly string TypeNameAttribute = "CLRTypeStored";
        public static readonly string CountAttribute = "CountStored";
        public static readonly string PropertiesAttribute = "PropertiesStored";

        // use Reflection.Emit to make methods that read and write the generic type T
        // This should be faster & allow us to remove the class constraint.
        // Currrently the class constraint exists because FieldInfo.SetValue doesn't work with value types.
        public static void WriteRecords<T>(IEnumerable<T> data, IGroup group) 
            //we remove this constraint to give callers more flexibility in usage
            //though we will do checks internally to see if we can process the type
            //where T : class, new()
        {
            data = data.ToList();
            var count = data.Count();

            // Add logic to control which fields go to HDF5 here
            var writeType = typeof(T);
            var props = writeType.GetProperties();

            var dataSpace = group.File.CreateDataspace(new long[] { count }, new long[] { count });

            IDataset[] datasets = new IDataset[props.Length];


            Func<PropertyInfo, int, Array> PrepFieldP = delegate(PropertyInfo info, int i)
            {
                var dataType = group.File.CreateDatatype(info.PropertyType);
                group.RemoveChild(info.Name);
                datasets[i] = group.CreateDataset(info.Name, dataType, dataSpace);
                return Array.CreateInstance(info.PropertyType, count);
            };

            Array[] dataArrays = props.Select(PrepFieldP).ToArray();

            // Use Reflection.Emit to make functions for copying data here
            int ec = 0;
            foreach (T el in data)
            {
                for (int j = 0; j < dataArrays.Length; j++)
                {
                    dataArrays[j].SetValue(props[j].GetValue(el, null), ec);
                }
                ec++;
            }

            // Transfer the C# arrays to the HDF5 file
            for (int j = 0; j < props.Length; j++)
            {
                datasets[j].Write(dataArrays[j]);
            }


            // Store some metadata in the IGroup
            // Store the name of the type T
            group.InsertAttribute(TypeNameAttribute, writeType.FullName);

            // Store the number of records written
            group.InsertAttribute(CountAttribute, count);

            // Store the list of field names stored
            group.InsertAttribute(PropertiesAttribute, String.Join(",", props.Select(f => f.Name).ToArray()));
        }

        public static IEnumerable<T> ReadRecords<T>(IGroup group)
        //we remove this constraint to give callers more flexibility in usage
        //though we will do checks internally to see if we can process the type
        //where T : class, new()
        {
            // Make sure the type we're loading matches the type stored in the group.
            var readType = typeof(T);
            var storedType = (string)group.GetAttribute(TypeNameAttribute).Read();

            if (storedType != readType.FullName)
            {
                throw new ArgumentException(String.Format("Attempted to read type '{0}' but found type '{1}' in IGroup '{2}' in Chunk File '{3}'",
                        readType.FullName, storedType, group.Name, group.File.Name));
            }

            var count = (int)group.GetAttribute(CountAttribute).Read();

            // The array of records that we'll return
            var resultArray = new T[count];

            for (int i = 0; i < count; i++)
            {
                T tObj = (T)Activator.CreateInstance(typeof (T));
                resultArray[i] = tObj;
            }

            // Get the set of fields that we need to copy into the structures.
            var storedFields = ((string)group.GetAttribute(PropertiesAttribute).Read()).Split(',');
            var copyFields = readType.GetProperties().Where(f => storedFields.Contains(f.Name)).ToArray();


            foreach (PropertyInfo f in copyFields)
            {
                Array data = (Array)((IDataset)group.GetChild(f.Name)).Read();
                for (int i = 0; i < count; i++)
                    f.SetValue(resultArray[i], data.GetValue(i), null);

            }

            return resultArray;
        }
    }


    /// <summary>
    /// Infrastructure for storing an array data classes to HDF5.  The user configures mappings
    /// from the storage type T to simple variables or array of simple variable that get written to 
    /// HDF5. 
    /// </summary>
    /// <typeparam name="T">The type to be stored in HDF5</typeparam>
    public class ArrayRecordStore<T>
    {
        public static int DatasetChunkSize = 16384;

        // Chunking for the ZMW and ZMWMetrics datasets.
        public static int ZMWChunkSizeX = 1024;
        public static int ZMWChunkSizeY = 256;

        // private static ILogText log = DiagManager.LogManager.LocalLogger();
        public static readonly string TypeNameAttribute = "CLRTypeStored";
        public static readonly string CountAttribute = "CountStored";
        public static readonly string PropertiesAttribute = "PropertiesStored";
        public static readonly string IndexDatasetName = "_Index_";
        public static readonly string DescriptionAttribute = "Description";
        public static readonly string UnitsAttribute = "UnitsOrEncoding";

        public static long[] LA(int i)
        {
            return new long[1] { i };
        }

        public static string GetIndexName(string name)
        {
            return String.Format("{0}.Index", name);
        }

        public delegate TQ[] ReaderFunc<TQ>(int start, int n, TQ[] inData);

        /// <summary>
        /// The writer class for persisting arrays of data with user-defined type T to HDF5.  You should 
        /// subclass this ArrayRecordStore&lt;T&gt;.Writer for the data type you are writing and call 
        /// the WriteSetupSingleton and WriteSetupArray methods to define 
        /// how the writer should pull data out of your datatype and store it in HDF5.
        /// </summary>
        public class Writer : IDisposable
        {
            private List<Action<T[]>> writers = new List<Action<T[]>>();
            private IGroup group;
            private int countStored = 0;
            internal uint compressionLevel = 0;
            internal bool compressData = false;

            public Writer(IGroup group, uint compressionLevel = 0u)
            {
                this.group = group;
                this.group.InsertAttribute(CountAttribute, countStored);

                if(compressionLevel > 0)
                {
                    this.compressionLevel = compressionLevel;
                    compressData = true;
                }
            }

            public IGroup CreateChildGroup(string name)
            {
                return group.CreateGroup(name);
            }

            /// <summary>
            /// The subclassing writer should call this method once for each simple field that needs to 
            /// be extracted from the data and stored in HDF5.
            /// </summary>
            /// <typeparam name="TR">The type of the field to be stored</typeparam>
            /// <param name="name">The name of the target HDF5 dataset</param>
            /// <param name="getter">A function that pulls out the relevant field from the data item being stored</param>
            /// <param name="outGroup">The group to store the dataset. Pass null to use the main Writer group</param>
            /// <param name="description">Optional descriptive text for the field</param>
            /// <param name="units"></param>
            public IDataset WriteSetupSingleton<TR>(string name, Func<T, TR> getter, IGroup outGroup = null,
                                                    string description = null, object units = null)
            {
                Type t = typeof(TR);
                var dataType = group.File.CreateDatatype(t);

                outGroup = outGroup ?? group;
                outGroup.RemoveChild(name);
                var dataSpace = outGroup.File.CreateDataspace(new long[] { 0 }, new long[] { -1 });

                var createPlist = new HDFDatasetCreateProperty();
                createPlist.Chunking = new long[] { ZMWChunkSizeX };

                // Turn on compression in the dataset if it is enabled.
                if (this.compressData)
                    createPlist.SetDeflate(this.compressionLevel);

                var dataset = ((HDFGroup)outGroup).CreateDataset(name, dataType, dataSpace, createPlist);
                
                //var dataset = outGroup.CreateDataset(name, dataType, dataSpace);

                if (description != null)
                {
                    dataset.InsertAttribute(DescriptionAttribute, description);
                }

                if (units != null)
                {
                    dataset.InsertAttribute(UnitsAttribute, units);
                }

                var size = 0;

                Action<T[]> writer = delegate(T[] data)
                {
                    var n = data.Length;
                    dataset.Extend(LA(size + n));
                    var ds = dataset.Dataspace;
                    ds.SelectHyperslab(LA(size), LA(1), LA(n), null);

                    TR[] outdata = data.Map(getter);
                    dataset.Write(outdata, ds);
                    size += n;
                };
                writers.Add(writer);

                return dataset;
            }


            /// <summary>
            /// The subclassing writer should call this method once for each simple field that needs to 
            /// be extracted from the data and stored in HDF5.
            /// </summary>
            /// <typeparam name="TR">The type of the field to be stored</typeparam>
            /// <param name="name">The name of the target HDF5 dataset</param>
            /// <param name="getter">A function that pulls out the relevant field from the data item being stored</param>
            /// <param name="outGroup">The group to store the dataset. Pass null to use the main Writer group</param>
            /// <param name="description">Optional descriptive text for the field</param>
            /// <param name="units"></param>
            public void WriteSetupFixedWidth<TR>(string name, Func<T, TR[]> getter, IGroup outGroup = null,
                                                 string description = null, object units = null)
            {
                Type t = typeof(TR);

                outGroup = outGroup ?? group;
                outGroup.RemoveChild(name);
                var dataSpace = outGroup.File.CreateDataspace(new long[] { 0, 1 }, new long[] { -1, -1 });

                // So, having allocated this IDisposable dataSet object and referenced it inside a
                // delegate that gets called who knows when or how many times... when does the dataSet
                // ever get disposed? Never?

                var createPlist = new HDFDatasetCreateProperty();
                createPlist.Chunking = new long[] { ZMWChunkSizeX, ZMWChunkSizeY };

                // Turn on compression in the dataset if it is enabled.
                if (this.compressData)
                    createPlist.SetDeflate(this.compressionLevel);

                IDataset dataSet;

                using (var dataType = group.File.CreateDatatype(t))
                    //dataSet = outGroup.CreateDataset(name, dataType, dataSpace);
                    dataSet = ((HDFGroup)outGroup).CreateDataset(name, dataType, dataSpace, createPlist);

                if (description != null)
                    dataSet.InsertAttribute(DescriptionAttribute, description);

                if (units != null)
                    dataSet.InsertAttribute(UnitsAttribute, units);

                var size = 0;
                int width = 0;

                Action<T[]> writer = delegate(T[] data)
                {
                    // This closure has been de-LINQ ified to make Mono happy.
                    // Apparently mono has some bugs in the compilation of code with numerous layers of 
                    // generics and closures.  I haven't pinned down the exact source though...

                    if (data.Length == 0)
                        return;

                    width = getter(data[0]).Length;
                    var outdata = new TR[data.Length, width];
                    var n = data.Length;

                    for (int i = 0; i < data.Length; i++)
                    {
                        var row = getter(data[i]);

                        if (row.Length != width)
                        {
                            throw new ArgumentException("Data array are not all the same width");
                        }

                        for (int j = 0; j < row.Length; j++)
                        {
                            outdata[i, j] = row[j];
                        }
                    }

                    dataSet.Extend(new long[] { size + n, width });
                    var ds = dataSet.Dataspace;

                    ds.SelectHyperslab(
                        new long[] { size, 0 },
                        new long[] { 1, 1 },
                        new long[] { n, width },
                        null);

                    dataSet.Write(outdata, ds);
                    size += n;
                };
                writers.Add(writer);
            }

            public class ArrayGroup
            {
                private readonly Writer writer;
                private readonly IGroup group;
                private readonly string indexFieldName;
                private readonly Dictionary<string, Func<T, int>> lengthCheckers = new Dictionary<string, Func<T, int>>();
                private readonly Dictionary<string, Action<int, T[]>> writers = new Dictionary<string, Action<int, T[]>>();
                private IDataset indexDataset;

                internal ArrayGroup(string indexName, Writer parent)
                {
                    writer = parent;
                    indexFieldName = indexName;
                    group = writer.group;
                }

                public void AddArrayField<TR>(string name, Func<T, IList<TR>> getter, IGroup outGroup =null,
                                              string description = null, object units = null)
                {
                    outGroup = outGroup ?? group;

                    Type t = typeof(TR);
                    var dataType = outGroup.File.CreateDatatype(t);

                    var hdfGroup = (HDFGroup) outGroup;
                    var createPlist = new HDFDatasetCreateProperty();
                    createPlist.Chunking = new long[] { DatasetChunkSize };

                    // Turn on compression in the dataset if it is enabled.
                    if(writer.compressData)
                        createPlist.SetDeflate(writer.compressionLevel);

                    // Create dataspaces to hold the raw data and the index array
                    outGroup.RemoveChild(name);
                    var dataSpaceData = hdfGroup.File.CreateDataspace(new long[] { 0 }, new long[] { -1 });
                    var dataDS = hdfGroup.CreateDataset(name, dataType, dataSpaceData, createPlist);

                    // Description annotation
                    if (description != null)
                        dataDS.InsertAttribute(DescriptionAttribute, description);

                    // Units annotation
                    if (units != null)
                        dataDS.InsertAttribute(UnitsAttribute, units);

                    // Annotate the new dataSpace with the name of the index field
                    dataDS.InsertAttribute("IndexField", indexFieldName);

                    lengthCheckers[name] = data => getter(data).Count;

                    Action<int, T[]> writerDelegate = delegate(int numExistingPoints, T[] data)
                    {
                        // See how many new things we're adding for each record
                        var fieldData = data.Map(getter);
                        var pointsPerRecord = fieldData.Map(a => a.Count);
                        var numNewPoints = pointsPerRecord.Sum();

                        // Make the array of new indices, and the flattened raw data.
                        var outdata = SafeStack<TR>(fieldData);

                        // Extend & write to the raw data dataset
                        dataDS.Extend(LA(numExistingPoints + numNewPoints));
                        var ds = dataDS.Dataspace;
                        ds.SelectHyperslab(LA(numExistingPoints), LA(1), LA(numNewPoints), null);
                        dataDS.Write(outdata, ds);
                    };

                    writers[name] = writerDelegate;
                    dataType.Dispose();     // Coverity RESOURCE_LEAK
                }

                
                public void AddArrayField<TR>(string name, Func<T, TR[,]> getter, IGroup outGroup = null,
                                              string description = null, object units = null)
                {
                    outGroup = outGroup ?? group;

                    Type t = typeof(TR);
                    var dataType = outGroup.File.CreateDatatype(t);

                    var hdfGroup = (HDFGroup) outGroup;

                    // Create dataspaces to hold the raw data and the index array
                    //var indexFieldName = getIndexName(name);
                    outGroup.RemoveChild(name);

                    // Setup chunking to keep the all the channels in a single chunk.
                    var createPlist = new HDFDatasetCreateProperty();
                    createPlist.Chunking = new long[] { DatasetChunkSize, 4 };

                    // Turn on compression in the dataset if it is enabled.
                    if (writer.compressData)
                        createPlist.SetDeflate(writer.compressionLevel);

                    var dataSpaceData = hdfGroup.File.CreateDataspace(new long[] { 0, 1 }, new long[] { -1, 4 });
                    var dataDS = hdfGroup.CreateDataset(name, dataType, dataSpaceData, createPlist);

                    // Description annotation
                    if (description != null)
                        dataDS.InsertAttribute(DescriptionAttribute, description);

                    // Units annotation
                    if (units != null)
                        dataDS.InsertAttribute(UnitsAttribute, units);

                    // Annotate the new dataSpace with the name of the index field
                    dataDS.InsertAttribute("IndexField", indexFieldName);

                    lengthCheckers[name] = data => getter(data).GetLength(0);

                    bool widthInit = false;
                    int width = 0;

                    Action<int, T[]> writerDelegate = delegate(int numExistingPoints, T[] data)
                    {
                        // See how many new things we're adding for each record
                        var fieldData = data.Map(getter);

                        if (widthInit == false)
                        {
                            width = fieldData[0].GetLength(1);
                            widthInit = true;
                        }

                        if (!fieldData.All(d => d.GetLength(1) == width))
                        {
                            throw new ArgumentException("Data array are not all the same width");
                        }


                        var pointsPerRecord = fieldData.Map(a => a.GetLength(0));
                        var numNewPoints = pointsPerRecord.Sum();

                        // Make the array of new indices, and the flattened raw data.
                        var outdata = SafeStack2(fieldData);

                        Func<int, long[]> lw = i => new long[] { i, width };

                        // Extend & write to the raw data dataset
                        dataDS.Extend(lw(numExistingPoints + numNewPoints));
                        var ds = dataDS.Dataspace;
                        ds.SelectHyperslab(
                            new long[] { numExistingPoints, 0 },
                            new long[] { 1, 1 },
                            new long[] { numNewPoints, width }, null);
                        dataDS.Write(outdata, ds);
                    };

                    writers[name] = writerDelegate;
                }
                
                /// <summary>
                /// Concatenate a series of arrays into a single array
                /// </summary>
                /// <param name="arrays">The input arrays</param>
                /// <returns>The concatenated arrays</returns>
                public static TA[] SafeStack<TA>(IList<TA>[] arrays)
                {
                    var l = arrays.Select(a => a.Count).Sum();

                    var result = SafeAlloc<TA>(l);
                    int count = 0;

                    foreach (var a in arrays)
                    {
                        for (int i = 0; i < a.Count; i++)
                            result[count + i] = a[i];

                        count += a.Count;
                    }

                    return result;
                }

                public static TA[,] SafeStack2<TA>(IList<TA[,]> data)
                {
                    if (data.Count == 0)
                        return new TA[0, 0];

                    var nRows = data.Select(i => i.GetLength(0)).Sum();
                    var nCols = data.Select(i => i.GetLength(1)).Max();

                    var res = SafeAlloc2<TA>(nRows, nCols);

                    int startRow = 0;
                    foreach (var d in data)
                    {
                        d.CopyTo(res, startRow, 0);
                        startRow += d.GetLength(0);
                    }

                    return res;
                }

                private static TR[,] SafeAlloc2<TR>(int rows, int cols)
                {
                    try
                    {
                        return new TR[rows,cols];
                    }
                    catch(OutOfMemoryException)
                    {
                        GC.Collect();
                        GC.WaitForPendingFinalizers();
                    }
                    // If this allocation fails, it will bring down the pipeline.  Too bad :(
                    return new TR[rows, cols];
                }

                private static TR[] SafeAlloc<TR>(int n)
                {
                    try
                    {
                        return new TR[n];
                    }
                    catch (OutOfMemoryException)
                    {
                        GC.Collect();
                        GC.WaitForPendingFinalizers();
                    }
                    // If this allocation fails, it will bring down the pipeline.  Too bad :(
                    return new TR[n];
                }


                private Action<T[]> BuildWriter(IGroup indexGroup)
                {
                    indexGroup.RemoveChild(indexFieldName);

                    using (var createPlist = new HDFDatasetCreateProperty())
                    {
                        createPlist.Chunking = new long[] {DatasetChunkSize};

                        using (var dataSpaceIndex = indexGroup.File.CreateDataspace(new long[] {0}, new long[] {-1}))
                        {
                            indexDataset = ((HDFGroup) indexGroup).CreateDataset(indexFieldName,
                                                                    indexGroup.File.CreateDatatype(typeof (int)),
                                                                    dataSpaceIndex, createPlist);

                            // Annotate for consistency w/ other markup fields
                            indexDataset.InsertAttribute(DescriptionAttribute, "ZMW event-stream counts");
                        }
                    }
                    // Keep track of how many things we've added
                    var numExistingRecords = 0;
                    var numExistingPoints = 0;

                    Action<T[]> groupWriter = data =>
                    {
                        // First check that all the fields have the same length for all the data items.
                        Func<T, bool> checkLengths = t =>
                        {
                            int n0 = lengthCheckers.First().Value(t);

                            foreach (var checker in lengthCheckers)
                            {
                                if (checker.Value(t) != n0)
                                    return false;
                            }

                            return true;
                        };

                        // Make sure each data item has consistent field lengths
                        if (!data.All(checkLengths))
                            throw new ArgumentException(
                                "Write failed because grouped fields don't have the same number of elements");

                        var getLength = lengthCheckers.First().Value;

                        var numNewRecords = data.Length;
                        var pointsPerRecord = data.Select(getLength).ToArray();
                        var numNewPoints = pointsPerRecord.Sum();

                        // Extend & write to the index dataset
                        indexDataset.Extend(
                            LA(numExistingRecords + numNewRecords));
                        var ids = indexDataset.Dataspace;
                        ids.SelectHyperslab(LA(numExistingRecords), LA(1),
                                            LA(numNewRecords), null);
                        indexDataset.Write(pointsPerRecord, ids);

                        // Write all the columns out
                        writers.Values.ForEach(w => w(numExistingPoints, data));

                        // Update the chunk pointers
                        numExistingRecords += numNewRecords;
                        numExistingPoints += numNewPoints;
                    };

                    return groupWriter;
                }

                public void Close()
                {
                    writer.writers.Add(BuildWriter(group));
                }

                public void Close(IGroup indexGroup)
                {
                    writer.writers.Add(BuildWriter(indexGroup));
                }
            }


            public ArrayGroup MakeParallelArrayGroup(string indexName)
            {
                return new ArrayGroup(indexName, this);
            }

            /// <summary>
            /// The subclassing writer should call this method for members of the storage datatype that 
            /// is a list of a primitive datatype.
            /// </summary>
            /// <typeparam name="TR">The primitive datatype to be stored for this field.</typeparam>
            /// <param name="name">The name of the target HDF5 dataset.</param>
            /// <param name="getter">A function that pulls out the relevant data array from the data item being stored</param>
            public void WriteSetupArray<TR>(string name, Func<T, IList<TR>> getter)
            {
                var ag = new ArrayGroup(GetIndexName(name), this);
                ag.AddArrayField(name, getter);
                ag.Close();
            }

            /// <summary>
            /// Write a series of data objects to the HDF5 file
            /// </summary>
            /// <param name="data">The data to be written</param>
            public void WriteRecords(IEnumerable<T> data)
            {
                WriteRecords(data, 1);
            }

            /// <summary>
            /// Write a series of data objects to the HDF5 file.  The data is broken into chunks for writing.
            /// </summary>
            /// <param name="data">The data to be written</param>
            /// <param name="chunkSize">The number of objects to write simultaneously in a single block</param>
            public void WriteRecords(IEnumerable<T> data, int chunkSize)
            {
                foreach (var chunk in data.Chunk(chunkSize))
                {
                    var chunkArray = chunk.ToArray();
                    writers.ForEach(w => w(chunkArray));
                    countStored += chunkArray.Length;
                    group.InsertAttribute(CountAttribute, countStored);
                }
            }

            public void AddExtraWriter(Action<T[]> w)
            {
                writers.Add(w);
            }

            #region IDisposable Members

            public void Dispose()
            {
                Dispose(true);
                GC.SuppressFinalize(this);
            }

            protected void Dispose(bool disposing)
            {
                if (disposing && group != null)
                {
                    group.Dispose();
                    group = null;
                }
            }

            ~Writer()
            {
                Dispose(false);
            }

            #endregion
        }

        public class Reader
        {
            private IGroup group;

            protected Reader(IGroup group)
            {
                this.group = group;
            }

            // Create a reader function for singleton properties
            protected SingletonReader<R> MakeSingletonReader<R>(string name)
            {
                // Type t = typeof(R);
                var dataset = (IDataset)group.GetChild(name);
                return new SingletonReader<R>(dataset);
            }

            // Create a reader function for indexed
            protected ArrayReader<R> MakeArrayReader<R>(string name)
            {
                // Type t = typeof(R);

                var dDataset = (IDataset)group.GetChild(name);

                string indexFieldName;
                IDataContainer idxField = dDataset.GetAttribute("IndexField");


                if (idxField != null)
                {
                    indexFieldName = (string)dDataset.GetAttribute("IndexField").Read();
                }
                else
                {
                    indexFieldName = GetIndexName(name);
                }


                var iDataset = (IDataset)group.GetChild(indexFieldName);

                // var numRecords = iDataset.Dataspace.Dimensions[0];
                var recordSizes = (int[])iDataset.Read();

                return new ArrayReader<R>(dDataset, recordSizes);
            }

            public int NumRecords { get { return (int)group.GetAttribute(CountAttribute).Read(); } }
        }
    }
}
