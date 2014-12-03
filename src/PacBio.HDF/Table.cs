using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Linq.Expressions;
using PacBio.Utils;

namespace PacBio.HDF
{
    public interface IStreamObserver<T> : IDisposable
    {
        void WriteStream(IEnumerable<T> data);
        void OnData(T[] data);
        void OnCompleted();
    }

    /// <summary>
    /// A class for storing simple C# types to parallel HDF5 arrays.
    /// </summary>
    public class Table<T> where T : new()
    {
        private static PacBioLogger log = PacBioLogger.GetLogger(2);

        static readonly string TypeNameAttribute = "CLRTypeStored";
        static readonly string CountAttribute = "CountStored";
        static readonly string PropertiesAttribute = "PropertiesStored";

        public static IStreamObserver<T> GetTableSink(string filename)
        {
            int bufferSize = 10000;
            return new Writer(filename, bufferSize);
        }


        public class Writer : IDisposable, IStreamObserver<T>
        {
            private SimpleTypeWriter[] writers;
            private IChunkFile f;

            public Writer(string filename, int bufferSize)
            {
                // Add logic to control which fields go to HDF5 here
                var writeType = typeof(T);
                var props = writeType.GetFields(BindingFlags.Instance | BindingFlags.Public | BindingFlags.NonPublic);

                f = HDFFile.Open(filename, FileMode.Create, FileAccess.ReadWrite);
                f.InsertAttribute(TypeNameAttribute, writeType.Name);

                writers = props.Map(p => new SimpleTypeWriter(p, bufferSize, p.Name, (HDFGroupContainer)f));

            }

            private int recordsWritten = 0;

            public void WriteChunk(T[] chunk)
            {
                var cl = chunk;
                recordsWritten += cl.Length;

                foreach (var w in writers)
                {
                    w.WriteChunk(cl);
                }

                f.InsertAttribute(CountAttribute, recordsWritten);
            }

            public void OnData(T[] data)
            {
                WriteChunk(data);
            }

            public void WriteStream(IEnumerable<T> data)
            {
                const int chunkSize = 2 << 16;
                var chunkStream = data.Chunk(chunkSize).Select(c => c.ToArray());
                var n = 0;

                foreach (var chunk in chunkStream)
                {
                    n += chunk.Length;
                    WriteChunk(chunk);
                    Console.WriteLine("Written {0}", n);
                }
            }

            public void OnCompleted()
            {
                // Close down the file
                Dispose();
            }

            protected virtual void Dispose(bool disposing)
            {
                if (disposing && f != null)
                {
                    f.Dispose();
                    f = null;
                }
            }

            public void Dispose()
            {
                Dispose(true);
                GC.SuppressFinalize(this);
            }

            ~Writer()
            {
                Dispose(false);
            }
        }


        /// <summary>
        /// Write data elements to HDF file
        /// </summary>
        /// <param name="filename">Target filename</param>
        /// <param name="data">Dataset to write</param>
        public static void WriteTable(string filename, IEnumerable<T> data)
        {
            int bufferSize = 10000;
            using (var w = new Writer(filename, bufferSize))
            {
                foreach (var chunk in data.Chunk(bufferSize))
                {
                    w.WriteChunk(chunk.ToArray());
                }
            }
        }

        /// <summary>
        /// Write data elements to HDF file
        /// </summary>
        /// <param name="filename">Target filename</param>
        /// <param name="data">Dataset to write</param>
        public static void WriteTable(string filename, IEnumerable<T[]> data)
        {
            int bufferSize = 10000;
            using (var w = new Writer(filename, bufferSize))
            {
                foreach (var chunk in data)
                {
                    w.WriteChunk(chunk);
                }
            }
        }

        public static IEnumerable<T[]> ReadTable(string filename)
        {
            return ReadTable(filename, 10000);
        }

        /// <summary>
        /// Read data a stream of data from an HDF5 table
        /// </summary>
        public static IEnumerable<T[]> ReadTable(string filename, int bufferSize)
        {
            // Add logic to control which fields go to HDF5 here
            var writeType = typeof(T);
            var props = writeType.GetFields(BindingFlags.Instance | BindingFlags.Public | BindingFlags.NonPublic);

            using (var f = HDFFile.Open(filename, FileMode.Open, FileAccess.Read))
            {
                if (f.GetAttribute(TypeNameAttribute) != null)
                {
                    var typeName = f.GetAttribute(TypeNameAttribute).ReadSingleton<string>();
                    if (typeName != writeType.Name)
                        throw new Exception(String.Format("Attempting to read type: {0}, but file contains type: {1}",
                                                          writeType.Name, typeName));
                }

                var readers = props.Map(p => new SimpleTypeReader(p, bufferSize, p.Name, (HDFGroupContainer)f));

                var total = f.GetAttribute(CountAttribute).ReadSingleton<int>();

                for (int i = 0; i < total; i += bufferSize)
                {
                    // Set up a chunk of the output type
                    var chunkSize = Math.Min(bufferSize, total - i);
                    var buf = new T[chunkSize];

                    for (int j = 0; j < chunkSize; j++)
                    {
                        buf[j] = new T();
                    }

                    foreach (var reader in readers)
                    {
                        reader.ReadChunk(buf, i);
                    }

                    yield return buf;

                }
            }
        }


        static long[] LA(int i)
        {
            return new long[1] { i };
        }


        private class SimpleTypeWriter
        {
            private Array buf;
            private IDataset dataset;

            private Func<int, Array, T[], int> reader;

            private int numExistingPoints = 0;
            private FieldInfo propInfo;
            private int bufferSize;

            public SimpleTypeWriter(FieldInfo p, int bufferSize, string dsName, HDFGroupContainer group)
            {
                // Setup chunking to keep the all the channels in a single chunk.
                var createPlist = new HDFDatasetCreateProperty();
                createPlist.Chunking = new long[] { bufferSize };

                // Turn on compression in the dataset
                createPlist.SetDeflate(3);

                this.bufferSize = bufferSize;
                propInfo = p;

                buf = Array.CreateInstance(p.FieldType, bufferSize);

                var dataSpace = group.File.CreateDataspace(new long[] { 0 }, new long[] { -1 });

                var dataType = group.File.CreateDatatype(p.FieldType);
                group.RemoveChild(dsName);
                dataset = group.CreateDataset(dsName, dataType, dataSpace, createPlist);

                reader = MakeReader(p);
            }

            public void WriteChunk(T[] data)
            {
                // Copy the new data from the raw chunk into the buffer
                var numNewPoints = data.Length;

                if (buf == null || numNewPoints != buf.Length)
                    buf = Array.CreateInstance(propInfo.FieldType, numNewPoints);

                reader(numNewPoints, buf, data);

                // Extend & write to the raw data dataset
                dataset.Extend(LA(numExistingPoints + numNewPoints));
                var ds = dataset.Dataspace;
                ds.SelectHyperslab(LA(numExistingPoints), LA(1), LA(numNewPoints), null);
                dataset.Write(buf, ds);

                numExistingPoints += numNewPoints;
            }

            public Func<int, Array, T[], int> MakeReader(FieldInfo p)
            {
                Func<int, Array, T[], int> f =
                    (n, bb, src) =>
                    {
                        for (int i = 0; i < n; i++)
                        {
                            object v = p.GetValue(src[i]);
                            bb.SetValue(v, i);
                        }
                        return n;
                    };

                return f;
            }

        }

        private class SimpleTypeReader
        {
            private Array buf;
            private IDataset dataset;

            private Func<Array, T[], int> reader;

            private FieldInfo propInfo;
            private int bufferSize;

            public SimpleTypeReader(FieldInfo p, int bufferSize, string dsName, HDFGroupContainer group)
            {
                this.bufferSize = bufferSize;
                propInfo = p;

                buf = Array.CreateInstance(p.FieldType, bufferSize);
                dataset = (IDataset)group.GetChild(dsName);
                reader = MakeReader(p);

            }

            public void ReadChunk(T[] data, int start)
            {
                // Copy the new data from the raw chunk into the buffer
                var numPoints = data.Length;

                if (buf == null || numPoints != buf.Length)
                    buf = Array.CreateInstance(propInfo.FieldType, numPoints);

                // Read the raw dataset
                var ds = dataset.Dataspace;
                ds.SelectHyperslab(LA(start), LA(1), LA(numPoints), null);
                dataset.Read(ref buf, ds);

                reader(buf, data);

            }

            private static Func<int, Array, T[], int> MakeReaderSlow(FieldInfo p)
            {
                Func<int, Array, T[], int> f =
                    (n, bb, target) =>
                    {
                        for (int i = 0; i < n; i++)
                        {
                            var val = bb.GetValue(i);
                            object tt = target[i];

                            p.SetValue(tt, val);
                            target[i] = (T)tt;
                        }

                        return n;
                    };

                return f;
            }

            public static Func<Array, T[], int> MakeReader(FieldInfo p)
            {
                var src = Expression.Parameter(typeof(Array), "src");
                var destArr = Expression.Parameter(typeof(T[]), "dest");

                var n = Expression.Variable(typeof(int), "n");

                var srcArrayType = p.FieldType.MakeArrayType();
                var srcGenericArr = Expression.Variable(srcArrayType, "srcGenericArr");

                var srcIdxExpr = Expression.ArrayIndex(srcGenericArr, n);
                var destArrIdxExpr = Expression.ArrayIndex(destArr, n);

                var assign = Expression.Assign(Expression.Field(destArrIdxExpr, p), srcIdxExpr);
                var incrExpr = Expression.PostIncrementAssign(n);
                var assignAndIncr = Expression.Block(assign, incrExpr);

                var label = Expression.Label(typeof(int));

                var loop =
                    Expression.Loop(
                        Expression.IfThenElse(Expression.LessThan(n, Expression.ArrayLength(srcGenericArr)),
                                              assignAndIncr,
                                              Expression.Break(label, n)),
                        label);

                var block = Expression.Block(
                    new[] { n, srcGenericArr },
                    Expression.Assign(n, Expression.Constant(0)),
                    Expression.Assign(srcGenericArr, Expression.TypeAs(src, srcArrayType)),
                    loop);

                Expression<Func<Array, T[], int>> l1 = Expression.Lambda<Func<Array, T[], int>>(
                    block, src, destArr);

                return l1.Compile();
            }
        }

        /// <summary>
        /// Write data records in the given HDF5 group 
        /// </summary>
        public static void WriteRecords(IEnumerable<T> data, IGroup group)
        // FIXME -- use Reflection.Emit to make methods that read and write the generic type T
        // This should be faster & allow us to remove the class constraint.
        // Currrently the class constraint exists because FieldInfo.SetValue doesn't work with value types.
        //we remove this constraint to give callers more flexibility in usage
        //though we will do checks internally to see if we can process the type
        //where T : class, new()
        {
            data = data.ToList();
            var count = data.Count();

            // Add logic to control which fields go to HDF5 here
            var writeType = typeof(T);
            var props = writeType.GetProperties();

            log.Log(LogLevel.DEBUG,
                    String.Format("Storing {0} instances of type '{1}' in IGroup '{2}'",
                                  count, writeType.FullName, group.Name));

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

        /// <summary>
        /// Read records
        /// </summary>
        public static IEnumerable<T> ReadRecords(IGroup group)
        //we remove this constraint to give callers more flexibility in usage
        //though we will do checks internally to see if we can process the type
        //where T : class, new()
        {
            // Make sure the type we're loading matches the type stored in the group.
            var readType = typeof(T);
            var storedType = (string)group.GetAttribute(TypeNameAttribute).Read();

            if (storedType != readType.FullName)
            {
                log.Log(LogLevel.ERROR,
                        String.Format("Attempted to read type '{0}' but found type '{1}' in IGroup '{2}' in Chunk File '{3}'",
                                      readType.FullName, storedType, group.Name, group.File.Name));
                throw new ArgumentException(String.Format("HDF Stored type '{0}' does not match request type '{1}",
                                                          storedType, readType.FullName));
            }

            var count = (int)group.GetAttribute(CountAttribute).Read();

            // The array of records that we'll return
            var resultArray = new T[count];

            for (int i = 0; i < count; i++)
            {
                T tObj = (T)Activator.CreateInstance(typeof(T));
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
}
