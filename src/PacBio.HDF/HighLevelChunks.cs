using System;
using System.IO;

using System.Linq;

namespace PacBio.HDF
{
    public class HighLevelChunks : IHighLevelChunks, IDisposable
    {
        public string FileName { get { return file.FileName; }}

        /// <summary>
        /// Use this method to access HighLevelChunks if you may be opening the same resource multiple times in the process.
        /// </summary>
        public static HighLevelChunks Open(string filename, bool forWrite)
        {
            return new HighLevelChunks(filename, forWrite);
        }

        public HighLevelChunks(string filename, bool forWrite)
        {
            var mode = forWrite ? FileMode.OpenOrCreate : FileMode.Open;
            var access = forWrite ? FileAccess.ReadWrite : FileAccess.Read;

            file = HDFFile.Open(filename,  mode, access);
        }



        protected IChunkFile file;

        /// <summary>
        /// Work with a specified (already open) IChunkFile
        /// Used when this is contained within another writer
        /// </summary>
        /// <param name="iChunkFile"></param>
        public HighLevelChunks(IChunkFile iChunkFile)
        {
            file = iChunkFile;
        }


        private static bool didInit = false;


        /// <summary>
        /// Open a chunk file based on parsing its URI
        /// </summary>
        /// <param name="filename"></param>
        /// <param name="forWrite"></param>
        /// <param name="access"></param>
        /// <returns></returns>
        public static IChunkFile Open(string filename, FileMode forWrite, FileAccess access)
        {
            return HDFFile.Open(filename, forWrite, access);
        }

        #region IChunkReader Members

        public void Close()
        {
            file.Dispose();
        }

         
        public object ReadDataset(string datasetName)
        {
            IDataset ds = (IDataset) file.GetChild(datasetName);

            if(ds == null)
                return null;        // oops - nothing here by that name

            return ds.Read();
        }


        /// <summary>
        /// Create an attribute on a specified node
        /// </summary>
        /// <param name="nodePath">The HDF path to the node (which must already exist)</param>
        /// <param name="attrName">The name of the attribute</param>
        /// <param name="attrValue">The value of the attribute - attribute type will be detected automatically</param>
        public IDataContainer WriteAttribute(string nodePath, string attrName, object attrValue)
        {
            IAttributeTarget ds = (IAttributeTarget)file.GetChild(nodePath);
            if (ds == null)
                throw new IOException("Node not found: " + nodePath);

            IDataContainer attr = ds.GetAttribute(attrName);
            Array arr = null;
            if (attr == null)
            {
                arr = attrValue as Array;
                if (arr == null)
                {
                    attr =
                        ds.CreateAttribute(attrName, file.CreateDatatype(attrValue.GetType()),
                                           file.CreateDataspace(new long[] {1}, new long[] {1}));
                }
                else
                {
                    // create the appropriate array-based attribute
                    long[] dims = new long[arr.Rank];
                    for (int i = 0; i < arr.Rank; i++)
                        dims[i] = arr.GetUpperBound(i) + 1;

                    attr =
                        ds.CreateAttribute(
                            attrName, file.CreateDatatype(attrValue.GetType()),
                            file.CreateDataspace(dims, dims));
                }
            }

            if (attr != null)
                attr.Write(attrValue);

            return attr;
        }

        /// <summary>
        /// Read an attribute on a specified node
        /// </summary>
        /// <param name="nodePath">the path to the node</param>
        /// <param name="attrName">the name of the attribute</param>
        /// <returns>null if not found</returns>
        public object ReadAttribute(string nodePath, string attrName)
        {
            using(var ds = (IAttributeTarget)file.GetChild(nodePath))
            {
                if (ds == null)
                    throw new IOException("Node not found: " + nodePath);

                IDataContainer attr = ds.GetAttribute(attrName);
                if (attr != null)
                    return attr.Read();
            }
            return null;
        }

        /// <summary>
        /// Read an attribute on a specified node
        /// </summary>
        /// <param name="nodePath">the path to the node</param>
        /// <param name="attrName">the name of the attribute</param>
        /// <returns>null if not found</returns>
        public T ReadAttributeSingleton<T>(string nodePath, string attrName)
        {
            using (var ds = (IAttributeTarget)file.GetChild(nodePath))
            {
                if (ds == null)
                    throw new IOException("Node not found: " + nodePath);

                IDataContainer attr = ds.GetAttribute(attrName);
                if (attr != null)
                {
                    var res = attr.Read();

                    if (res is T)
                        return (T)res;

                    if (res is T[])
                        return ((T[])res)[0];
                }
            }

            throw new ArgumentException("Attribute not found");
        }


        #endregion

        #region IChunkWriter Members


        public void WriteDataset(string datasetName, object value)
        {
            ChunkUtils.Create(file, datasetName, value);
        }

        #endregion

        #region IDisposable Members

        protected virtual void Dispose(bool disposing)
        {
            if (disposing && file != null)
            {
                file.Dispose();
                file = null;
            }
        }

        public void Dispose()
        {
            Dispose(true);
            GC.SuppressFinalize(this);
        }

        ~HighLevelChunks()
        {
            Dispose(false);
        }

        #endregion

        #region IHighLevelChunks Members

        public IChunkFile File
        {
            get { return file; }
        }

        public IDataspace CreateDataspace(int[] dims, int[] maxDims)
        {
            return file.CreateDataspace(dims.Cast<long>().ToArray(), maxDims.Cast<long>().ToArray());
        }

        public IDataset WriteDataset(string datasetName, object value, IDataspace dspace)
        {
            // Does it already exist?
            IDataset dset = (IDataset) file.GetChild(datasetName);

            if (dset == null)
            {
                using (IDatatype dtype = file.CreateDatatype(value.GetType()))
                {
                    IGroup mygroup = ChunkUtils.ExpandPath(file, ref datasetName);

                    dset = mygroup.CreateDataset(datasetName, dtype, dspace);
                }
            }

            dset.Write(value, dspace);

            return dset;
        }



        public object ReadDataset(string datasetName, IDataspace dspace)
        {
            // If we don't have a dataspace specified, then just return the 'standard' Read result
            if (dspace == null)
                return ReadDataset(datasetName);

            IDataset ds = (IDataset)file.GetChild(datasetName);
            if(ds == null)
                return null;

            Array result = Array.CreateInstance(ds.Datatype.NativeType, dspace.NumElements);

            ds.Read(ref result, dspace);

            return result;
        }

        /// <summary>
        /// Create a group (or open an existing group)
        /// </summary>
        /// <param name="path"></param>
        /// <returns></returns>
        public IGroup CreateGroup(string path)
        {
            IChunkElement e = file.GetChild(path);

            if(e != null)
            {
                IGroup g = e as IGroup;

                if (g == null)
                    throw new IOException("HDF path already a dataset: " + path);

                return g;
            }

            IGroup mygroup = ChunkUtils.ExpandPath(file, ref path);

            return mygroup.CreateGroup(path);
        }

        #endregion
    }
}
