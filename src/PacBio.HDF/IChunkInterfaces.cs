using System;

namespace PacBio.HDF
{

    /// <summary>
    /// Represent an HDF5 file.
    /// </summary>
    public interface IChunkFile : IGroup
    {
        /// <summary>
        /// A factory to get a data type description given a native .net type
        /// </summary>
        /// <param name="t"></param>
        /// <returns></returns>
        IDatatype CreateDatatype(Type t);


        /// <summary>
        /// Create an array dataspace
        /// </summary>
        /// <param name="curDimensions">the size in each dimension</param>
        /// <param name="maxDimensions">the size in each dimension (or -1 for no maximum)</param>
        /// <returns></returns>
        IDataspace CreateDataspace(long[] curDimensions, long[] maxDimensions);


        /// <summary>
        /// Create a simple non array dataspace
        /// </summary>
        /// <returns></returns>
        IDataspace CreateDataspace();

        /// <summary>
        /// Uri of the HDF5 file
        /// </summary>
        string FileName { get; }

    }

    public interface IChunkElement : IDisposable
    {
        string Name { get; }

        /// <summary>
        /// Get the file that this element is part of
        /// </summary>
        IChunkFile File { get; }
    }


    /// <summary>
    /// Any object that can be tagged with attributes
    /// </summary>
    public interface IAttributeTarget : IChunkElement
    {
        /// <summary>
        /// Create an attribute
        /// </summary>
        /// <param name="name"></param>
        /// <param name="datatype"></param>
        /// <param name="dataspace"></param>
        /// <returns></returns>
        /// 
        IDataContainer CreateAttribute(string name, IDatatype datatype, IDataspace dataspace);

        /// <summary>
        /// Read a named attribute
        /// </summary>
        /// <param name="name"></param>
        /// <returns>null for not found</returns>
        IDataContainer GetAttribute(string name);

        /// <summary>
        /// Read all of the attributes attached to this object
        /// </summary>
        /// <returns></returns>
        IDataContainer[] GetAttributes();

        /// <summary>
        /// Deletes an attribute attached to this object
        /// </summary>
        /// <param name="name">Name of the attribute to delete</param>
        /// <returns>True if an attribute was found and deleted, false if not</returns>
        bool DeleteAttribute(string name);
    }


    /// <summary>
    /// Behavior common to anything that can actually contain real data 
    /// </summary>
    /// In the case of HDF that means Attribute or Datasets.
    public interface IDataContainer : IChunkElement
    {
        IDatatype Datatype { get; }
        IDataspace Dataspace { get; }

        /// <summary>
        /// Read the entire contents of this data chunk
        /// </summary>
        /// <returns></returns>
        object Read();

        /// <summary>
        /// Read a subsection of this chunk by using hyperslabs
        /// </summary>
        /// <param name="target">The array to fill with data from the chunk.  Read() may or may not replace the target object you pass in</param>
        /// <param name="dataspace">a dataspace defining what hyperslab we want to read, or null to read the whole hyperslab from the file</param>
        void Read(ref Array target, IDataspace dataspace);

        /// <summary>
        /// Write the entire contents of this chunk
        /// </summary>
        /// <param name="o"></param>
        void Write(object o);

        /// <summary>
        /// Write a subsection of this chunk
        /// </summary>
        /// <param name="o"></param>
        /// <param name="dataspace"></param>
        /// Use when you've created a dataspace with a restricted hyperslab
        void Write(object o, IDataspace dataspace);
    }


    public interface IDatatype : IChunkElement
    {
        /// <summary>
        /// The .net type that is appropriate for storing this datatype
        /// </summary>
        Type NativeType { get; }
    }

    /// <summary>
    /// Provides links into HDF5 "Region Reference" - particular hyperslab sections of datasets
    /// </summary>  
    public interface IReference : IChunkElement
    {
        /// <summary>
        /// Get the dataset that this reference is pointing to
        /// </summary>
        /// <returns></returns>
        IDataset Dereference();


        /// <summary>
        /// Get the selected region in the target dataset 
        /// </summary>
        /// <returns></returns>
        IDataspace GetRegion();
    }


    public interface IDataspace : IDisposable
    {
        /// <summary>
        /// return the # of elements in each dimension
        /// </summary>
        /// <returns>An array of length equal to the rank of this dataset, each element shows the size in that dimension</returns>
        long[] Dimensions { get; }

        /// <summary>
        /// The maximum amount a specified dimension may ever reach
        /// </summary>
        /// Or -1 for no limit
        long[] MaxDimensions { get; }


        /// <summary>
        /// Select a subportion of this dataspace
        /// </summary>
        /// <param name="start">the initial coordinates to use for the top corner of the hyperslab</param>
        /// <param name="stride">how much to increment each dim by as we move between blocks, or null for a stride of all 1s</param>
        /// <param name="count">the number of blocks in each direction</param>
        /// <param name="block">the dimensions of each block we are selecting, or null for blocks of size 1</param>
        /// If you call this method, the dataset it is attached to will use 'partial' IO for reads and writes.  Useful for
        /// stuff like movie frames.
        void SelectHyperslab(
            long[] start,
            long[] stride,
            long[] count,
            long[] block);

        /// <summary>
        /// Undo the effects of SelectHyperslab - go back to reading/writing the whole space
        /// </summary>
        void SelectAll();

        /// <summary>
        /// The # of dimensions for each block written to a hyperslab (or Dimensions if not using hyperslabs)
        /// </summary>
        long [] HyperBlock
        { 
            get;
        }
                /// <summary>
        /// How many elements in this dataset?
        /// </summary>
        /// Takes into account any hyperslabbing going on
        long NumElements
        {
            get;
        }
    }


    public interface IDataset : IDataContainer, IAttributeTarget
    {
        /// <summary>
        /// Expand a resizable dataset
        /// </summary>
        /// <param name="newDims"></param>
        void Extend(long[] newDims);
        
        /// <summary>
        /// Create a 'region reference' that points inside this particular dataset
        /// </summary>
        /// <param name="dataspace"></param>
        /// <returns></returns>
        /// It is the callers responsiblity to write this to the file via a dataset
        IReference CreateReference(IDataspace dataspace);
    }

    /// <summary>
    /// Any object that can contain groups or datasets
    /// </summary>
    public interface IGroup : IAttributeTarget
    {

        /// <summary>
        /// Create (or reopen an existing) group
        /// </summary>
        /// <param name="groupName"></param>
        /// <returns></returns>
        IGroup CreateGroup(string groupName);

        /// <summary>
        /// Create a dataset with a more flexible representation
        /// </summary>
        /// <param name="datasetName"></param>
        /// <param name="dataspace"></param>
        /// <param name="datatype"></param>
        IDataset CreateDataset(string datasetName, IDatatype datatype, IDataspace dataspace);

        /// <summary>
        /// Open an existing dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>null for not found</returns>
        IChunkElement GetChild(string datasetName);

        /// <summary>
        /// Remove a child dataset or group from a group
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>false if no child was deleted</returns>
        bool RemoveChild(string datasetName);

        /// <summary>
        /// Return the children of this node
        /// </summary>
        /// <returns></returns>
        IChunkElement[] GetChildren();
    }
}
