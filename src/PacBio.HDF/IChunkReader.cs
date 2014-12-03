using System;

namespace PacBio.HDF
{
    /// <summary>
    /// Operations for reading chunks
    /// </summary>
    public interface IChunkReader : IDisposable
    {
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
