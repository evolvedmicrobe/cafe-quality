using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Runtime.InteropServices;

namespace PacBio.HDF
{
    /// <summary>
    /// A high level API for reading and writing from any chunk provider (HDF, XML, AVI, ChunkPath etc...)
    /// </summary>
    /// The regular IChunkElement tree is designed to be a very small API to make it easy to implement new providers.
    /// This interface (and its associated implementation HighLevelChunks) provides an easier API for USERS of the chunk
    /// framework.  It makes it easy to insert/extract data without having to worry about datatypes, datasets etc...
    /// I recommend this as your primary gateway for opening chunk files
    public interface IHighLevelChunks
    {
        /// <summary>
        /// Return the low level file for advanced users
        /// </summary>
        /// May not be necessary
        IChunkFile File { get; }

        /// <summary>
        /// Return the Uri used to open the IHighLevelChunks object.
        /// </summary>
        string FileName { get; }

        /// <summary>
        /// Create an attribute on a specified node
        /// </summary>
        /// <param name="nodePath">The HDF path to the node (which must already exist)</param>
        /// <param name="attrName">The name of the attribute</param>
        /// <param name="attrValue">The value of the attribute - attribute type will be detected automatically</param>
        IDataContainer WriteAttribute(string nodePath, string attrName, object attrValue);

        /// <summary>
        /// Read an attribute on a specified node
        /// </summary>
        /// <param name="nodePath">the path to the node</param>
        /// <param name="attrName">the name of the attribute</param>
        /// <returns>null if not found</returns>
        object ReadAttribute(string nodePath, string attrName);

        IDataspace CreateDataspace( int[] dims, int[] maxDims);
        IDataset WriteDataset(string datasetName, object value, IDataspace dspace);
      
        object ReadDataset(string datasetName, IDataspace dspace);


        /// <summary>
        /// Create a group (or open an existing group)
        /// </summary>
        /// <param name="path"></param>
        /// <returns></returns>
        IGroup CreateGroup(string path);

        /// <summary>
        /// store a chunk in the file - if the named dataset does not exist it will be created
        /// </summary>
        /// <param name="datasetName">A four character type code</param>
        /// <param name="value"></param>
        /// Note: for legacy purposes we currently map some four character 'avi style' dataset names to new 
        /// HDF style names internally.  kevinh will add more details on this process eventually.
        /// This simple method assumes a simple 1 element dataset of either string, numeric or array type
        void WriteDataset(string datasetName, object value);


        /// <summary>
        /// Read the entire named dataset - returning data in a natural format
        /// </summary>
        /// <param name="datasetName">The dataset name used to fetch the specified dataset</param>
        /// <returns>the data type depends on how the data was defined at creation, or null if not found</returns>
        /// Note: for legacy purposes we currently map some four character 'avi style' dataset names to new 
        /// HDF style names internally.  kevinh will add more details on this process eventually.
        // [return: MarshalAs(UnmanagedType.SafeArray)]
        object ReadDataset(string datasetName);

    }
}
