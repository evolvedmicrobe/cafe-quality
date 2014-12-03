using System;

namespace PacBio.HDF
{
    /// <summary>
    /// Operations for writing chunks
    /// </summary>
    public interface IChunkWriter : IDisposable
    {
        /// <summary>
        /// store a chunk in the file - if the named dataset does not exist it will be created
        /// </summary>
        /// <param name="datasetName">A four character type code</param>
        /// <param name="value"></param>
        /// Note: for legacy purposes we currently map some four character 'avi style' dataset names to new 
        /// HDF style names internally.  kevinh will add more details on this process eventually.
        /// This simple method assumes a simple 1 element dataset of either string, numeric or array type
        void WriteDataset(string datasetName, object value);
    }


}
