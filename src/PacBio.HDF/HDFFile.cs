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
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text;
using System.Runtime.InteropServices;
using System.Text.RegularExpressions;
using PacBio.Utils;

#pragma warning disable 168

namespace PacBio.HDF
{
    public abstract class HDFChunkElement : IChunkElement
    {
        /// <summary>
        /// Free a block of memory alloced by the HDF library
        /// </summary>
        /// <param name="ptr"></param>
        protected static void HDFFree(IntPtr ptr)
        {
            // Use the cruntime free call - because the hdf5 routine used malloc
            CRuntime.free(ptr);
        }

        #region IChunkElement Members

        public abstract string Name 
        {
            get; 
        }

        #endregion


        public override string ToString()
        {
            return Name;
        }

        #region IDisposable Members

        protected abstract void Dispose(bool disposing);

        public void Dispose()
        {
            Dispose(true);
            GC.SuppressFinalize(this);
        }

        ~HDFChunkElement()
        {
            Dispose(false);
        }

        #endregion

        /// <summary>
        /// Our HDF5 identifier
        /// </summary>
        internal int Id { get; set; }


        /// <summary>
        /// Sometimes we might not have a valid Id yet, in that case throw an exception
        /// </summary>
        protected void CheckValid()
        {
            if (Id == 0)
                throw new ApplicationException("ID is invalid");
        }


        public abstract IChunkFile File
        { 
            get;
        }
    }


    /// <summary>
    /// Any HDF object that is inside of a file (i.e. not including HDFFile itself)
    /// </summary>
    public abstract class HDFIntraFile : HDFChunkElement, IAttributeTarget
    {
        private IChunkFile file;

        internal HDFIntraFile(IChunkFile parent)
        {
            file = parent;    
        }

        public IDataContainer GetAttribute(string name)
        {
            if(!HDFAttribute.Exists(this, name))
                return null;

            return new HDFAttribute(this, name);
        }

        public virtual IDataContainer[] GetAttributes()
        {
            int numattr = H5A.get_num_attrs(Id);

            var result = new IDataContainer[numattr];
            for (uint i = 0; i < numattr; i++)
                result[i] = new HDFAttribute(this, i);

            return result;
        }

        public IDataContainer CreateAttribute(string name, IDatatype datatype, IDataspace dataspace)
        {
            return new HDFAttribute(this, name, datatype, dataspace);
        }

        public bool DeleteAttribute(string name)
        {
            if (H5A.exists(Id, name))
            {
                H5A.delete(Id, name);
                return true;
            }
            return false;
        }

        public override string Name
        {
            get
            {
                CheckValid(); 
                
                return H5I.get_name(Id);
            }
        }

        public override IChunkFile File
        {
            get { return file; }
        }

        public int Copy(IGroup destinationGroup, string destinationName)
        {
            var hdfGroup = destinationGroup as HDFGroupContainer;

            if (hdfGroup == null)
                throw new Exception("Must copy HDF item to HDFGroup");

            return H5O.copy(Id, Name, hdfGroup.Id, destinationName);
        }

        public int Copy(IGroup destinationGroup)
        {
            return Copy(destinationGroup, Name);
        }
    }
    
    /// <summary>
    /// Encapsulate an H5G HDF group, and provides access to its children
    /// </summary>
    public class HDFGroup : HDFGroupContainer, IGroup
    {
        internal HDFGroup(HDFGroupContainer parent, string name) : base(parent.File)
        {
            Id = H5G.create(parent.Id, name);
        }

        /// <summary>
        /// Open an existing group 
        /// </summary>
        /// <param name="file"></param>
        /// <param name="id"></param>
        /// assuming the caller has already opened the specified id, and we are now responsible for closing it
        internal HDFGroup(IChunkFile file, int id) : base(file)
        {
            Id = id;
        }

        private bool disposed = false;

        protected override void Dispose(bool disposing)
        {
            if (!disposed)
            {
                disposed = true;
                H5G.close(Id);
            }
        }
    }

    /// <summary>
    /// Encapsulates the H5T datatype system, and provide mappings between .NET types and HDF5 types.
    /// </summary>
    public class HDFDatatype : HDFChunkElement, IDatatype
    {
        internal HDFDatatype(int id)
        {
            Id = id;
        }

        internal HDFDatatype(Type t)
        {
            Id = TypeToID(t);
        }

        static HDFDatatype()
        {
            InitTypemap(typeof (byte), "H5T_STD_U8LE");
            InitTypemap(typeof(sbyte), "H5T_NATIVE_SCHAR");
            InitTypemap(typeof(ulong), "H5T_NATIVE_ULLONG");
            InitTypemap(typeof(long), "H5T_NATIVE_LLONG");
            InitTypemap(typeof(ushort), "H5T_NATIVE_USHORT");
            InitTypemap(typeof(short), "H5T_NATIVE_SHORT");
            InitTypemap(typeof(uint), "H5T_NATIVE_UINT");
            InitTypemap(typeof(int), "H5T_NATIVE_INT");
            InitTypemap(typeof(float), "H5T_NATIVE_FLOAT"); 
            InitTypemap(typeof(double), "H5T_NATIVE_DOUBLE");

	    InitTypemap(typeof(H5R.RegRef), RefString);

            InitTypemap(typeof (string),
                        "H5T_STRING { STRSIZE H5T_VARIABLE ; STRPAD H5T_STR_NULLTERM ; CSET H5T_CSET_ASCII ; CTYPE H5T_C_S1; }");
        }

        private static int H5T_STRING = 3;
        private static int H5T_COMPOUND = 6;

        /// <summary>
        /// We special case the type lookups for this type, because it is busted in the standard HDF library
        /// </summary>
        private static int RefDtype;

        private const string RefString = "H5T_STD_REF_DSETREG";

        static void InitTypemap(Type t, string hdfname)
        {
            int dtype;

            // This type is not listed in the high level text_to_dtype function
            if (hdfname == RefString)
                dtype = RefDtype = HDFGlue.GetHDFGlobal(hdfname + "_g");
            else 
                dtype = H5L.text_to_dtype(hdfname);

            typeToHDF.Add(t, dtype);

            // We always search using the cononical string encoding
            string nativetype = DtypeToText(dtype);
            hdfToType.Add(nativetype, t);
        }

        static string DtypeToText(int dtype)
        {
            if(RefDtype != 0 && H5T.equal(dtype, RefDtype))
                return RefString;
            else
                return H5L.dtype_to_text(dtype);
        }

        public override string ToString()
        {
            return HDFTypestring;
        }

        /// <summary>
        /// The human readable cononical encoding for this datatype
        /// </summary>
        internal string HDFTypestring
        {
            get
            {
                CheckValid();

                return DtypeToText(Id);
            }
        }
       
        /// <summary>
        /// Type mappings for standard types
        /// </summary>
        private static IDictionary<Type, int> typeToHDF = new Dictionary<Type, int>();

        /// <summary>
        /// Type mappings for standard types
        /// </summary>
        /// Converts from an HDF native type to our .net types
        private static IDictionary<string, Type> hdfToType = new Dictionary<string, Type>();


        /// <summary>
        /// Convert a .net type to a HDF type code
        /// </summary>
        /// <param name="t"></param>
        /// <returns></returns>
        /// Creates a _new_ datatype object, the caller is responsible for freeing it
        static int TypeToID(Type t)
        {
            int hdfstr;

            // If someone is trying to write out our public HDFReference type, the data
            // stored in the file is actually the RegRef structure
            if (t == typeof(HDFReference))
                t = typeof (H5R.RegRef);    

            if (typeToHDF.TryGetValue(t, out hdfstr))
            {
                hdfstr = H5T.copy(hdfstr);

                return hdfstr;
            }

            throw new ApplicationException("Unknown native type to HDF: " + t);
        }

        /// <summary>
        /// The HDF type ID that is most efficient when storing this data type in RAM
        /// </summary>
        internal HDFDatatype MapToNativeType()
        {
            int nativetype = H5T.get_native_type(Id);
            return new HDFDatatype(nativetype);
        }

        /// <summary>
        /// A regex to find the length of a string
        /// </summary>
        static readonly Regex findLen = new Regex(@"STRSIZE\s+(?<digits>\d+)",
            RegexOptions.Compiled | RegexOptions.IgnoreCase);


        /// <summary>
        /// If this is a value type, how many bytes are needed to store an element?
        /// </summary>
        /// If this is a reference type, return 0.
        /// Strings are funny because HDF strings are passed by value if we know the length in advance, by ref otherwise
        public int ValueSize
        {
            get
            {
                if (NativeType == typeof(string))
                {
                    string tstring = HDFTypestring;

                    Match match = findLen.Match(tstring);
                    if (match.Success)
                    {
                        string strlen = match.Result("${digits}");
                        int blen;

                        if (int.TryParse(strlen, out blen))
                            return blen;
                    }

                    // Must be a variable length string - fall in and handle like any other reference type
                }

                return !NativeType.IsValueType ? 0 : Marshal.SizeOf(NativeType);
            }
        }

        private Type nativeCache = null;

        /// <summary>
        /// The .net type that is appropriate for storing this datatype
        /// </summary>
        public Type NativeType
        {
            get
            {
                if (nativeCache == null)
                {
                    foreach(var vp in typeToHDF)
                    {
                        if (H5T.equal(vp.Value, Id))
                            nativeCache = vp.Key;
                    }

                    if (nativeCache == null)
                    {
                        int hdfClass = H5T.get_class(Id);

                        if (hdfClass == H5T_STRING)
                            nativeCache = typeof (string);
                    }

                    if (nativeCache == null)
                    {
                        // We always search using native types  
                        Type ntype;
                        string tstring = HDFTypestring;

                        // Handle strings as a special format - because HDF has a zillion different encodings
                        if (tstring.StartsWith("H5T_STRING"))
                            nativeCache = typeof (string);
                        else if (hdfToType.TryGetValue(tstring, out ntype))
                            nativeCache = ntype;
                        else
                            throw new ApplicationException("Unknown HDF to native: " + HDFTypestring);
                    }
                }

                return nativeCache;
            }
        }

        public override string Name
        {
            get { return "DataType"; }
        }

        /// <summary>
        /// How many bytes are needed to hold this type?
        /// </summary>
        public long NumBytes
        {
            get
            {
                CheckValid();

                return H5T.get_size(Id);
            }
        }

        protected override void Dispose(bool disposing)
        {
            if (Id > 0)
                H5T.close(Id);

            Id = 0;
        }

        public override IChunkFile File
        {
            get { throw new NotImplementedException(); }
        }
    }

 
    public class HDFReference : HDFChunkElement, IReference
    {
        /// <summary>
        /// Our actual reference bytes
        /// </summary>
        private H5R.RegRef refptr;

        private HDFFile file;

        internal HDFReference(IDataset _target, IDataspace _dataspace)
        {
            HDFDataset target = (HDFDataset) _target;
            HDFDataspace dataspace = (HDFDataspace) _dataspace;

            file = (HDFFile) target.File;
            refptr = H5R.create(target.Id, ".", dataspace.Id);
        }

        /// <summary>
        /// Private constructor for upconverting from raw RegRefs
        /// </summary>
        /// <param name="file"></param>
        /// <param name="refptr"></param>
        internal HDFReference(HDFFile file, H5R.RegRef refptr)
        {
            this.file = file;
            this.refptr = refptr;
        }

        /// <summary>
        /// Our actual reference bytes
        /// </summary>
        internal H5R.RegRef RefPtr { get { return refptr; } }

        #region IReference Members

        public IDataset Dereference()
        {
            int newid = H5R.dereference(file.Id, RefPtr);

            return new HDFDataset(File, newid);
        }

        public IDataspace GetRegion()
        {
            int newid = H5R.get_region(file.Id, RefPtr);

            return new HDFDataspace(File, newid);
        }

        #endregion

        public override string Name
        {
            get { return "(HDFReference)"; }
        }

        protected override void Dispose(bool disposing)
        {
            // References live inside of datasets, so they have no resources that need to be freed.
        }

        public override IChunkFile File
        {
            get { return file; }
        }
    }


    public class HDFDataspace : HDFChunkElement, IDataspace
    {
        private IChunkFile parent;

        /// <summary>
        /// Create an array dataspace
        /// </summary>
        /// <param name="parent"></param>
        /// <param name="curDimensions">the size in each dimension</param>
        /// <param name="maxDimensions">the size in each dimension (or -1 for no maximum or null for same as dimensions)</param>
        /// <returns></returns>
        internal HDFDataspace(IChunkFile parent, long[] curDimensions, long[] maxDimensions) 
        {
            this.parent = parent;

            //This allow zero size of extent if maxDimension is H5S_UNLIMITED
            Id = H5S.create_simple(curDimensions, maxDimensions);  
        }



        /// <summary>
        /// Create a simple non array dataspace
        /// </summary>
        /// <param name="parent"></param>
        /// <returns></returns>
        internal HDFDataspace(IChunkFile parent) 
        {
            this.parent = parent;

            Id = H5S.create(H5S.SCALAR);
        }

        /// <summary>
        /// Create a dspace object given a dspace id
        /// </summary>
        /// <param name="id"></param>
        /// <param name="parent"></param>
        internal HDFDataspace(IChunkFile parent, int id)
        {
            this.parent = parent;

            Id = id;
        }

#if false
        /// <summary>
        /// Given a .net type return a dataspace that describes the in memory representation of that type
        /// </summary>
        /// <param name="t"></param>
        /// <param name="parent"></param>
        internal HDFDataspace(IChunkFile parent, Type t)
        {
            this.parent = parent;

            Id = H5S.copy(H5S.ALL);
        }
#endif


        public override string ToString()
        {
            var sb = new StringBuilder();

            sb.Append(Name);

            if (Dimensions.Length != 0)
            {
                sb.Append(" dim=[");

                foreach (long l in Dimensions)
                    sb.AppendFormat("{0} ", l);

                sb.Append("]");
                sb.Append(" numelem=" + NumElements);
            }

            return sb.ToString();
        }


        private long[] block = null;

        /// <summary>
        /// The # of dimensions for each block written to a hyperslab (or Dimensions if not using hyperslabs)
        /// </summary>
        public long [] HyperBlock
        {
            get { return block ?? Dimensions; }
        }


        /// <summary>
        /// The total size needed by the hyperslab (used to grow dimensions in dataset write.
        /// </summary>
        private long[] extent = null;

        /// <summary>
        /// The total size needed by the hyperslab (used to grow dimensions in dataset write)
        /// </summary>
        public long [] Extent
        {
            get { return extent ?? Dimensions; }
        }

        /// <summary>
        /// Select a subportion of this dataspace
        /// </summary>
        /// <param name="start">the initial coordinates to use for the top corner of the hyperslab</param>
        /// <param name="stride">how much to increment each dim by as we move between blocks, or null for a stride of all 1s</param>
        /// <param name="count">the number of blocks in each direction</param>
        /// <param name="block">the dimensions of each block we are selecting, or null for blocks of size 1</param>
        public void SelectHyperslab(
            long[] start,
            long[] stride,
            long[] count,
            long[] block)
        {
            HDFFile f = (HDFFile) File;

            this.block = block;

            if (block == null)
                block = start.Length.Fill(i => 1L);

            if (stride == null)
                stride = start.Length.Fill(i => 1L);

            extent = new long[start.Length];
            for (int i = 0; i < start.Length; i++)
                // extent = start + (count - 1) * stride + block;
                extent[i] = start[i] + (count[i] - 1) * stride[i] + block[i];

            H5S.select_hyperslab(Id, start, stride, count, block);
        }

        /// <summary>
        /// Undo the effects of SelectHyperslab - go back to reading/writing the whole space
        /// </summary>
        public void SelectAll()
        {
            block = null;
            extent = null;

            H5S.select_all(Id);
        }

        /// <summary>
        /// How many elements in this dataset?
        /// </summary>
        /// Takes into account any hyperslabbing going on
        public long NumElements
        {
            get
            {
                CheckValid();

                return H5S.get_select_npoints(Id);
            }
        }
        public override string Name
        {
            get { return (Dimensions.Length == 0) ? "ScalarDS" : "SimpleDS"; }
        }

        protected override void Dispose(bool disposing)
        {
            if (Id > 0)
                H5S.close(Id);

            Id = 0;
        }

        #region IDataspace Members

        long[] curdim, maxdim;

        void ReadDims()
        {
            H5S.get_simple_extent(Id, out curdim, out maxdim);
        }



        public long[] Dimensions
        {
            get
            {
                CheckValid();

                ReadDims();

                return curdim;
            }

            set
            {
                CheckValid();

                ReadDims();
                curdim = value;

                H5S.set_extent_simple(Id, curdim, maxdim);
            }
        }

        public long[] MaxDimensions
        {
            get
            {
                CheckValid();

                ReadDims();

                return maxdim;
            }
        }

        #endregion

        public override IChunkFile File
        {
            get { return parent; }
        }
    }


    public abstract class HDFDataContainer : HDFIntraFile, IDataContainer
    {
        internal HDFDataContainer(IChunkFile parent) : base(parent)
        {
        }


        #region IDataContainer Members

 
        public abstract IDatatype Datatype
        { 
            get;
        }

        public abstract IDataspace Dataspace
        { 
            get;
        }

        /// <summary>
        /// Create an array with the correct # of dimensions and datatype to read our data into
        /// </summary>
        /// <returns></returns>
        /// <param name="hdfRaw">If true, we want an array suitable for the HDF low level read/writes, if false we want an array suitable for C# users</param>
        /// Different hdfRaw settings will change this from expecting raw arrays of bytes or pointers to bytes
        Array CreateArray(bool hdfRaw)
        {
            // If the user asked for 'string' results we temporariy create an array appropriate for the type of HDF
            // read we are going to do
            int valSize;
            Type nativeType;

            using (var hdfType = (HDFDatatype)Datatype)
            {
                nativeType = hdfType.NativeType;
                valSize = hdfType.ValueSize;
            }

            // Assume we don't need to scale element sizes
            int elemSize = 1;
            
            long[] numdims;
            
            using (var space = Dataspace)
                numdims = space.Dimensions;

            // For string or other reference types we might need to substite an array of ptrs or an array of bytes
            if (hdfRaw && !nativeType.IsPrimitive)  // was !nativeType.IsValueType
            {
                if (valSize == 0)
                    nativeType = typeof (IntPtr);
                else
                {
                    // kevinh - we now support structs
                    // Debug.Assert(nativeType == typeof (string));        // We only understand strings now
                    nativeType = typeof (byte);
                    elemSize = valSize;
                    numdims = new long[] { elemSize * numdims.Aggregate(1L, (i, j) => i * j) };
                    //numdims[0] *= elemSize;
                }                
            }

            if (numdims.Length == 0)
            {
                // Must be a scalar dataspace, just claim we want a one element array
                numdims = new long[] { 1 };                
            }

            /*
            if(elemSize != 1 && hdfRaw)
            {
                Debug.Assert(numdims.Length == 1);      // We only understand one dim string arrays right now
                numdims[0] *= elemSize;
            }
            */

            Array arr = Array.CreateInstance(nativeType, numdims);
            return arr;
        }


        private readonly static UTF8Encoding fromAscii = new UTF8Encoding();


        /// <summary>
        /// Read a subsection of this chunk by using hyperslabs
        /// </summary>
        /// <param name="target">The array to fill with data from the chunk</param>
        /// <param name="filedataspace">a dataspace defining what data we want to read</param>
        /// Like Read but we assume that any array dimension swapping has already been done
        void ReadRowMajor(ref Array target, IDataspace filedataspace)
        {
            // The datatype we should use to match the layout of .net objects
            using (var datasetType = (HDFDatatype) Datatype)
            using (var hdfNativeType = datasetType.MapToNativeType())
            {
                // We try to have HDF do any necessary int size conversion for us,
                // but we can't do that if we are reading pointer based objects (arrays of strings etc)
                Type outputArrayType = target.GetType().GetElementType();

                // We can not automatically convert _to IntPtrs_ or _from strings/refs_.
                Type fromType = hdfNativeType.NativeType;
                bool isSimpleType = (outputArrayType != typeof(IntPtr)) && (fromType != typeof(string)) && (fromType != typeof(H5R.RegRef));

                // Size our memory dataspace by looking at the dimensions of the array that was passed in
                Tuple<Type, HDFDataspace> tu = GetElementType(target, filedataspace);
                using (var memdataspace = tu.Item2)
                using (HDFDatatype dtypeInRAM = isSimpleType ? new HDFDatatype(outputArrayType) : new HDFDatatype(H5T.copy(hdfNativeType.Id)))
                {
                    // Read our data into a 1d array of bytes
                    Array linear = target;

                    // Don't let this buffer move
                    GCHandle pinArray = GCHandle.Alloc(linear, GCHandleType.Pinned);

                    try
                    {
                        IntPtr bytes = Marshal.UnsafeAddrOfPinnedArrayElement(linear, 0);

                        ReadLow(dtypeInRAM, memdataspace, (HDFDataspace) filedataspace, bytes);

                        // Copy our bytes into our structured array
                        // No longer needed now that we read in place
                        // Buffer.BlockCopy(linear, 0, target, 0, nbytes);

                        // If our array is an array of pointers, we'll need to create objects of the correct type and copy those objects over
                        Type desiredType = dtypeInRAM.NativeType;

                        if (!desiredType.IsPrimitive)   // Was !IsValueType
                        {
                            Array newArray = CreateArray(false);

                            int srcoffset = 0;
                            int elemSize = dtypeInRAM.ValueSize;

                            foreach (long[] index in newArray.EnumerateIndicies()) // i = 0; i < newArray.Length; i++)
                            {
                                object newval;

                                if (desiredType == typeof(string))
                                {
                                    if (elemSize == 0)   // We must be using ptrs
                                    {
                                        var inarg = (IntPtr) target.GetValue(srcoffset);
                                        newval = Marshal.PtrToStringAnsi(inarg);
                                        HDFFree(inarg);
                                    }
                                    else
                                    {
                                        // copy from in place
                                        var btarget = (byte[]) target;
                                        // cut off the tring with the null terminator -- crap was getting appended
                                        var str = fromAscii.GetString(btarget, srcoffset, elemSize);

                                        if (str.IndexOf((char) 0) > 0)
                                            newval = str.Substring(0, str.IndexOf((char) 0));
                                        else
                                            newval = str;
                                    }
                                }
                                else
                                {
                                    Debug.Assert(desiredType.IsValueType);  // If not a string, the type better be a struct
                                    IntPtr srcptr = Marshal.UnsafeAddrOfPinnedArrayElement(linear, srcoffset);
                                    newval = Marshal.PtrToStructure(srcptr, desiredType);
                                }

                                newArray.SetValue(newval, index);
                                srcoffset += (elemSize != 0) ? elemSize : 1;
                            }
                            target = newArray;
                        }
                    }
                    finally
                    {
                        pinArray.Free();
                    }
                }
            }
        }

        /// <summary>
        /// Read a subsection of this chunk by using hyperslabs
        /// </summary>
        /// <param name="target">The array to fill with data from the chunk</param>
        /// <param name="filedataspace">a dataspace defining what data we want to read</param>
        public void Read(ref Array target, IDataspace filedataspace)
        {
#pragma warning disable 168
            using (var leakCheck = new HDFLeakChecker((HDFFile) File))
#pragma warning restore 168
            {
                if (filedataspace != null)
                {
                    // Grow our dataset dims to be at least as large as the file dimensions (the user will use the hyperslab
                    // for smaller reads)
                    using (var space = Dataspace)
                        ((HDFDataspace)filedataspace).Dimensions = space.Dimensions;
                }

                if (target.Length > 0)
                {
                    ReadRowMajor(ref target, filedataspace);
                }
            }
        }


        public object Read()
        {
            using (var leakCheck = new HDFLeakChecker((HDFFile)File))
            {
                Array array = CreateArray(true);

                ReadRowMajor(ref array, null);

                object result = array;

                // If the dataspace is actually for a scalar field, just return the first field of the array
                using (var dspace = Dataspace)
                    if (dspace.Dimensions.Length == 0)
                        result = array.GetValue(0);

                // If the result is a raw hdf reference ptr, convert it to our high level class
                if (result is H5R.RegRef)
                    result = new HDFReference((HDFFile)File, (H5R.RegRef)result);

                return result;
            }
        }

        public override string ToString()
        {
            using (var dt = Datatype)
            using (var ds = Dataspace)
                return base.ToString() + " " + dt + " " + ds;
        }

        /// <summary>
        /// The low level read operation
        /// </summary>
        /// <param name="memtype"></param>
        /// <param name="memdspace"></param>
        /// <param name="filedspace"></param>
        /// <param name="rawBytes">This object must already be allocated with enough empty bytes to contain the expected read</param>
        protected abstract void ReadLow(HDFDatatype memtype, HDFDataspace memdspace, HDFDataspace filedspace, IntPtr rawBytes);

        /// <summary>
        /// The low level write operation (slightly different glue for the dataset/attribute cases
        /// </summary>
        /// <param name="memtype"></param>
        /// <param name="memdspace"></param>
        /// <param name="filedspace"></param>
        /// <param name="buf"></param>
        protected abstract void WriteLow(HDFDatatype memtype, HDFDataspace memdspace, HDFDataspace filedspace, IntPtr buf);

        static UTF8Encoding utf8Converter = new UTF8Encoding();


        /// <summary>
        /// The MarshalObject call may need to allocate different types of memory, we group that data through this MemRef
        /// </summary>
        /// You must dispose this mem ref to ensure memory is released/unpinned
        class MemRef : IDisposable
        {
            public MemRef(IntPtr primary, GCHandle? primaryHandle, IntPtr secondary, MemoryPool<byte>.Ref poolRef)
            {
                this.primaryPtr = primary;
                this.handle = primaryHandle;
                this.secondary = secondary;
                this.poolRef = poolRef;
            }

            /// <summary>
            /// A low level ptr passed to the HDF lib
            /// </summary>
            IntPtr primaryPtr;

            /// <summary>
            /// If we need to alloc an extra array (for strings, we also need to free this)
            /// </summary>
            IntPtr secondary;

            /// <summary>
            /// Our primaryPtr might be a raw Ptr or have an associated GC handle that also needs to be disposed
            /// </summary>
            GCHandle? handle;

            /// <summary>
            /// If the memory came from a pool we'll also need to free this
            /// </summary>
            MemoryPool<byte>.Ref poolRef;

            /// <summary>
            /// The actual ptr needed by the HDF libraries
            /// </summary>
            public IntPtr Contents
            {
                get { return primaryPtr; }
            }

            #region IDisposable Members

            private bool disposed = false;

            protected virtual void Dispose(bool disposing)
            {
                if (!disposed)
                {
                    disposed = true;

                    if (secondary != IntPtr.Zero)
                        Marshal.FreeHGlobal(secondary);

                    // Super crufty - since UnsafeAddrOfPinnedArray element doesn't return a 'handle' we don't want to call free
                    // on ocopy if we have been provided with a gchandle to use instead
                    if (handle.HasValue)
                        handle.Value.Free();
                    else
                        Marshal.FreeHGlobal(primaryPtr);

                    poolRef = null; // Allow any pool memory to be GCed
                }
            }

            public void Dispose()
            {
                Dispose(true);
                GC.SuppressFinalize(this);
            }

            ~MemRef()
            {
                Dispose(false);
            }

            #endregion
        }

        /// <summary>
        /// Copy an object to memory that matches the HDF writer's expectations
        /// </summary>
        /// <param name="o"></param>
        /// <returns></returns>
        /// The caller must free the handle with FreeHGlobal!!!
        MemRef MarshalObject(object o)
        {
            // For references, we just write the RefPtr structure inside the object
            HDFReference refobj = o as HDFReference;
            if (refobj != null)
                o = refobj.RefPtr;

            var str = o as string;
            if (str != null)
            {
                // HDF _sometimes_ expects strings to be passed in as an _array of str ptrs_ so we need to add a second level of indirection
                // If the strings are FIXED length in the HDF definition, then it expects a single array of packed chars
                int strlen;

                using (HDFDatatype hdfType = (HDFDatatype)Datatype)
                {
                    strlen = hdfType.ValueSize;
                }

                if (strlen == 0) // Must be passing strings by reference
                {
                    var secondaryPtr = Marshal.StringToHGlobalAnsi(str);
                    int ptrsize = Marshal.SizeOf(secondaryPtr);
                    IntPtr buf = Marshal.AllocHGlobal(ptrsize);
                    Marshal.WriteIntPtr(buf, secondaryPtr);

                    return new MemRef(buf, null, secondaryPtr, null);
                }
                else
                {
                    // Get a ptr to the bytes
                    return new MemRef(Marshal.StringToHGlobalAnsi(str), null, IntPtr.Zero, null);
                }
            }
            else if(o is Array)
            {
                var arr = (Array) o;

                Type elementType = arr.GetType().GetElementType();

                if (elementType == typeof(string))
                {
                    // HDF expects us to pass in an array of ptrs for
                    // things like arrays of strings.  

                    int nelem = arr.Length;

                    var linear = new IntPtr[nelem];

                    // Preflight to find how many bytes we need to store the temporary copy of the strings (possibly allocating too many bytes)
                    int bytesNeeded = 0;
                    foreach (long[] indices in arr.EnumerateIndicies()) 
                    {
                        // Copy the element in the array
                        var element = (string)arr.GetValue(indices);

                        // We assume we need to add a 0 byte terminator after each str
                        bytesNeeded += utf8Converter.GetMaxByteCount(element.Length) + 1;
                    }

                    // We store our string contents in the secondary ptr (the ptrs to the strings go in the return value)
                    var secondaryPtr = Marshal.AllocHGlobal(bytesNeeded);

                    int ptrNum = 0;     // the ptr we are currently building
                    int charNum = 0;    // the byte # we are filling in our packed array of 0 terminated strs

                    // Hmm - We need to collapse a possibly 2d array to 1d array of char ptrs
                    // Collapse to 1-d array
                    foreach (long[] indices in arr.EnumerateIndicies()) 
                    {
                        // Copy the element in the array
                        var element = (string) arr.GetValue(indices);

                        byte[] bytes = utf8Converter.GetBytes(element);

                        // copy the str bytes
                        int curStrOffset = charNum;
                        Marshal.Copy(bytes, 0, secondaryPtr.Offset(charNum), bytes.Length);
                        charNum += bytes.Length;

                        // Add a 0 term
                        Marshal.WriteByte(secondaryPtr.Offset(charNum), 0);
                        charNum++;

                        linear[ptrNum++] = secondaryPtr.Offset(curStrOffset);
                    }

                    // Don't let this buffer move
                    var secondaryHandle = GCHandle.Alloc(linear, GCHandleType.Pinned);
                    return new MemRef(Marshal.UnsafeAddrOfPinnedArrayElement(linear, 0), secondaryHandle, secondaryPtr, null);
                }
                else
                {
                    // Don't let this buffer move
                    //var secondaryHandle = GCHandle.Alloc(linear, GCHandleType.Pinned);
                    var secondaryHandle = GCHandle.Alloc(arr, GCHandleType.Pinned);

                    //return new MemRef(Marshal.UnsafeAddrOfPinnedArrayElement(linear, 0), secondaryHandle, IntPtr.Zero, poolref);
                    return new MemRef(Marshal.UnsafeAddrOfPinnedArrayElement(arr, 0), secondaryHandle, IntPtr.Zero, null);
                }
            }
            else // Some other value type
            {
                int sz = Marshal.SizeOf(o);
                IntPtr buf = Marshal.AllocHGlobal(sz);
                Marshal.StructureToPtr(o, buf, false);

                return new MemRef(buf, null, IntPtr.Zero, null);
            }
        }

        /// <summary>
        /// Return the type of elements in the given array (or if not an array just return the object).
        /// This allows us to treat single elements like arrays of one element.
        /// </summary>
        Tuple<Type, HDFDataspace> GetElementType(object obj, IDataspace fileDataspace)
        {
            Type t = obj.GetType();

            if (t.IsArray)
            {
                Type nativeType;
                int elemSize;

                using (HDFDatatype hdftype = ((HDFDatatype)Datatype))
                {
                    nativeType = hdftype.NativeType;
                    elemSize = hdftype.ValueSize;
                }

                // It might be that we want to scale other array dims,
                // but for now I'm sure that struct arrays are messed up -
                // so be careful to only change those.

                bool isStruct = !nativeType.IsPrimitive && elemSize != 0;
                if (!isStruct)
                    elemSize = 1;
                
                var arr = (Array) obj;

                var dims = arr.Rank.Fill(i => arr.GetLongLength(i)/elemSize);

                // The native H5Dread and H5Dwrite methods do work when the dimensions 
                // (rank) of the memory and file dataspace differ (assuming that the number
                // of elements agree), but they run much more slowly.  
                // So coerce the dspace to agree in dimensionality
                // with the fileDataspace, not the target array.

                if (fileDataspace != null && fileDataspace.Dimensions.Length > arr.Rank)
                {
                    dims = GetOptimalMemspaceDims(fileDataspace,arr,elemSize);
                }
                    
                var maxDims = dims.Map(v => v > 0 ? v : -1);

                return new Tuple<Type, HDFDataspace>(t.GetElementType(), new HDFDataspace(File, dims, maxDims));
            }
            else
            {
                // must be a scalar

                return new Tuple<Type, HDFDataspace>(t, null);
            }
        }

        private static long[] GetOptimalMemspaceDims(IDataspace fileDataspace, Array arr, int elemSize)
        {
            int rankFile = fileDataspace.Dimensions.Length;
            int rankArr = arr.Rank;

            var dims = rankFile.Fill(1L);

            // Assign the dims consistent with the target array and rank of the
            // file dataspace by filling from the fastest-changing dimension in
            // the target array.  This may produce a memory dataspace with dimensions
            // that differ from the target hyperslab.  But, we assume that the
            // memory dataspace (in this use) is always tightly packed (non hyperslab),
            // and so ordering is preserved.
            //
            for (int i = 1; i <= rankArr; i++)
            {
                dims[rankFile - i] = arr.GetLength(rankArr - i)/elemSize;
            }
                       
            return dims;
        }

        public virtual void Write(object o, IDataspace filedataspace)
        {
            Tuple<Type, HDFDataspace> tu = GetElementType(o, filedataspace);

            using (var memdataspace = tu.Item2)
            using (HDFDatatype t = new HDFDatatype(tu.Item1))
            using (var ocopy = MarshalObject(o))
                WriteLow(t, memdataspace, (HDFDataspace) filedataspace, ocopy.Contents);
        }


        public void Write(object o)
        {
            Write(o, null);
        }

        /// <summary>
        /// The target will be filled with the data from filedataspace
        /// </summary>
        /// <returns></returns>
        
        public void CompoundTypeRead(HDFDatatype dtype, ref object target, IDataspace filedataspace)
        {
            Array arr = (Array)target;

            long[] dims = new long[1];
            dims[0] = arr.GetLength(0);
            
            int sz = Marshal.SizeOf(arr.GetType().GetElementType());
            int nbytes = arr.Length * sz;

            byte[] linear = new byte[nbytes];
            IntPtr ptr = Marshal.AllocHGlobal(sz);
            IntPtr bytes = Marshal.UnsafeAddrOfPinnedArrayElement(linear, 0);
            
            // Size our memory dataspace by looking at the dimensions of the array that was passed in
            using(var memdataspace = new HDFDataspace(File, dims, dims))
                ReadLow(dtype, memdataspace, (HDFDataspace)filedataspace, bytes);

            for (int i = 0; i < arr.Length; i++)
            {
                Marshal.Copy(linear, i * sz, ptr, sz);
                arr.SetValue(Marshal.PtrToStructure(ptr, arr.GetValue(i).GetType()), i);
            }

            Marshal.FreeHGlobal(ptr);
        }
     



        public void CompoundTypeWrite(HDFDatatype t, object o, IDataspace memdataspace, IDataspace filedataspace)
        {
            using (var ocopy = MarshalObject(o))
            {
                WriteLow(t, (HDFDataspace)memdataspace, (HDFDataspace)filedataspace, ocopy.Contents);
            }
        }

        public void CompoundTypeWrite(HDFDatatype t, object o, IDataspace filedataspace)
        {
            using (var ocopy = MarshalObject(o))
            {
                WriteLow(t, null, (HDFDataspace)filedataspace, ocopy.Contents);
            }
        }

        public void CompoundTypeWrite(HDFDatatype t, object o)
        {
            CompoundTypeWrite(t, o, null);
        }

 

        #endregion
    }


    public class HDFAttribute : HDFDataContainer, IDataContainer
    {
        /// <summary>
        /// Create a new attribute
        /// </summary>
        /// <param name="parent"></param>
        /// <param name="_name"></param>
        /// <param name="datatype"></param>
        /// <param name="dataspace"></param>
        internal HDFAttribute(IAttributeTarget parent, string _name, IDatatype datatype, IDataspace dataspace) : base(parent.File)
        {
            Id = H5A.create(((HDFChunkElement) parent).Id, _name, ((HDFDatatype)datatype).Id, ((HDFDataspace)dataspace).Id);
        }

        /// <summary>
        /// Open an existing attribute
        /// </summary>
        /// <param name="parent"></param>
        /// <param name="_name"></param>
        internal HDFAttribute(IAttributeTarget parent, string _name)
            : base(parent.File)
        {
            Id = H5A.open(((HDFChunkElement)parent).Id, _name);
        }


        /// <summary>
        /// Look up a child attribute by 0 based index
        /// </summary>
        /// <param name="parent"></param>
        /// <param name="index"></param>
        internal HDFAttribute(IAttributeTarget parent, uint index)
            : base(parent.File)
        {
            Id = H5A.open_idx(((HDFChunkElement) parent).Id, index);
        }


        /// <summary>
        /// Does the named attribute exist?
        /// </summary>
        /// <param name="parent"></param>
        /// <param name="_name"></param>
        /// <returns></returns>
        internal static bool Exists(IAttributeTarget parent, string _name)
        {
            return H5A.exists(((HDFChunkElement)parent).Id, _name);
        }
       

        public override IDatatype Datatype
        {
            get
            {
                CheckValid(); 
                
                return new HDFDatatype(H5A.get_type(Id));
            }
        }

        public override IDataspace Dataspace
        {
            get
            {
                CheckValid();

                return new HDFDataspace(File, H5A.get_space(Id));
            }
        }

        public override string Name
        {
            get
            {
                return H5A.get_name(Id);
            }
        }

        protected override void Dispose(bool disposing)
        {
            if (Id > 0)
                H5A.close(Id);

            Id = 0;
        }

        protected override void WriteLow(HDFDatatype memtype, HDFDataspace memdspace, HDFDataspace filedspace, IntPtr buf)
        {
            Debug.Assert(filedspace == null);       // subchunking aint supported

            H5A.write(Id, memtype.Id, buf);
        }

        protected override void ReadLow(HDFDatatype memtype, HDFDataspace memdspace, HDFDataspace filedspace, IntPtr rawBytes)
        {
            Debug.Assert(filedspace == null);       // subchunking aint supported

            H5A.read(Id, memtype.Id, rawBytes);
        }
    }



    public class HDFProperty : HDFChunkElement
    {
        private readonly string proptype;

        public HDFProperty(string proptype)
        {
            this.proptype = proptype;
            Id = H5P.create(proptype);
        }

        public HDFProperty(string proptype, int id)
        {
            this.proptype = proptype;
            Id = id;
        }

        public override string Name
        {
            get { return proptype; }
        }

        protected override void Dispose(bool disposing)
        {
            if (Id > 0)
                H5P.close(Id);

            Id = 0;
        }

        public override IChunkFile File
        {
            get { throw new NotImplementedException(); }
        }
    }

    public class HDFDatasetCreateProperty : HDFProperty
    {
        public HDFDatasetCreateProperty()
            : base("H5P_CLS_DATASET_CREATE")
        {
        }

        public HDFDatasetCreateProperty(int id) : base("H5P_CLS_DATASET_CREATE", id)
        {
        }

        /// <summary>
        /// The HDF chunking for datasets that can grow later
        /// </summary>
        public long [] Chunking
        {
            set
            {
                H5P.set_chunk(Id, value);
            }
        }

        public long[] GetChunking(int ndims)
        {
            return H5P.get_chunk(Id, ndims);
        }

        /// <summary>
        /// Turn on GZip compression for this dataset.  
        /// </summary>
        /// <param name="level">The gzip compression level, The level must be between 0-9 inclusive.  A lower level is
        /// faster, but gives less compression.</param>
        public void SetDeflate(uint level)
        {
            if(level < 0 || level > 9)
                throw new Exception("Bad GZIP compression level");


            H5P.set_deflate(Id, level);
        }
    }

    public class HDFFileAccessProperty : HDFProperty
    {
        public HDFFileAccessProperty() : base("H5P_CLS_FILE_ACCESS_ID")
        {
        }
        
        public void SetCache(int nMetaDataElements, int nChunksCached, ulong chunkCacheBytes, double preemptFactor)
        {
            H5P.set_cache(Id, nMetaDataElements, nChunksCached, chunkCacheBytes, preemptFactor);
        }
    }

    public class HDFDatasetAccessProperty : HDFProperty
    {
        public HDFDatasetAccessProperty()
            : base("H5P_CLS_DATASET_ACCESS")
        {
        }

        public void SetCache(uint nChunksCached, uint chunkCacheBytes, double preemptFactor)
        {
            H5P.set_chunk_cache(Id, nChunksCached, chunkCacheBytes, preemptFactor);
        }
    }
    


    public class HDFDataset : HDFDataContainer, IDataset
    {
        /// <summary>
        /// Create a new dataset
        /// </summary>
        /// <param name="parent"></param>
        /// <param name="datasetName"></param>
        /// <param name="dataspace"></param>
        /// <param name="datatype"></param>
        internal HDFDataset(HDFGroupContainer parent, string datasetName, IDatatype datatype, IDataspace dataspace)
            : base(parent.File)
        {
            // Needed to support variable length data sets
            using(HDFDatasetCreateProperty prop = new HDFDatasetCreateProperty())
            {
                if (!Extensions.ArrayEquals(dataspace.Dimensions, dataspace.MaxDimensions))
                {
                    long[] chunks = new long[dataspace.Dimensions.Length];

                    for (int i = 0; i < chunks.Length; i++)
                        // Assume the chunks are large for the dimensions that are not allowed to change
                        //chunks[i] = (dataspace.Dimensions[i] == dataspace.MaxDimensions[i]) ? dataspace.Dimensions[i] : 1;
                        //chunks on all dims (extendable or not) are set the same as Dimensions
                        if (dataspace.Dimensions[i] != 0)
                        {
                            chunks[i] = dataspace.Dimensions[i];
                        }
                        else
                        {
                            chunks[i] = 1024;
                        }

                    prop.Chunking = chunks;
                }
 
                Id = H5D.create(parent.Id, datasetName, ((HDFDatatype) datatype).Id, ((HDFDataspace) dataspace).Id,
                                H5P.DEFAULT, prop.Id);
            }
        }

        /// <summary>
        /// Create a new dataset, and allow the caller to specify properties
        /// </summary>
        /// <param name="parent"></param>
        /// <param name="datasetName"></param>
        /// <param name="datatype"></param>
        /// <param name="dataspace"></param>
        /// <param name="prop"></param>
        internal HDFDataset(HDFGroupContainer parent, string datasetName, IDatatype datatype,
                            IDataspace dataspace, HDFDatasetCreateProperty prop)
            : base(parent.File)
        {
                Id = H5D.create(parent.Id, datasetName, ((HDFDatatype)datatype).Id,
                                ((HDFDataspace)dataspace).Id, H5P.DEFAULT, prop.Id);
        }

        /// <summary>
        /// Open an existing dataset
        /// </summary>
        /// <param name="parent"></param>
        /// <param name="datasetName"></param>
        internal HDFDataset(HDFGroupContainer parent, string datasetName)
            : base(parent.File)
        {
            Id = H5D.open(parent.Id, datasetName);
        }


        /// <summary>
        /// Open an existing dataset
        /// </summary>
        /// <param name="parent"></param>
        /// <param name="datasetName"></param>
        /// <param name="plist"></param>
        internal HDFDataset(HDFGroupContainer parent, string datasetName, HDFDatasetAccessProperty plist)
            : base(parent.File)
        {
            Id = H5D.open(parent.Id, datasetName, plist.Id);
        }


        /// <summary>
        /// Open an existing dataset
        /// </summary>
        /// <param name="file"></param>
        /// <param name="id"></param>
        /// Note: we assume the ID has already been opened some other place and we are now responsible for closing it.
        internal HDFDataset(IChunkFile file, int id) : base(file)
        {
            Id = id;
        }


        /// <summary>
        /// Expand a resizable dataset
        /// </summary>
        /// <param name="newDims"></param>
        public void Extend(long[] newDims)
        {
            H5D.extend(Id, newDims);
        }


        public override void Write(object o, IDataspace filedataspace)
        {
            using (var leakCheck = new HDFLeakChecker((HDFFile)File))
            {
                // We automatically grow the dataset extents on disk if needed
                if (filedataspace != null)
                {
                    HDFDataspace fspace = (HDFDataspace)filedataspace;

                    Extend(fspace.Extent);

                    // Now that we grew our dataset, we might need to update the callers dataspace
                    using (var ds = Dataspace)
                    {
                        fspace.Dimensions = ds.Dimensions;
                    }
                }

                base.Write(o, filedataspace);
            }
        }

        protected override void Dispose(bool disposing)
        {
            if (Id > 0)
                H5D.close(Id);

            Id = 0;
        }

        #region IDataContainer Members

        /// <summary>
        /// A cached copy of the Datasets Datatype
        /// </summary>
        public override IDatatype Datatype
        {
            get 
            { 
                CheckValid();
                return new HDFDatatype(H5D.get_type(Id));
            }
        }

        /// <summary>
        /// A fresh copy of the Datasets Dataspace
        /// </summary>
        public override IDataspace Dataspace
        {
            get
            {
                CheckValid();
                return new HDFDataspace(File, H5D.get_space(Id));
            }
        }


        #endregion

        #region IAttributeTarget Members

        /* These should all b live in HDFIntraFile
        public IAttribute CreateAttribute(string name, IDatatype datatype, IDataspace dataspace)
        {
            return new HDFAttribute(this, name, datatype, dataspace);
        }

        public bool DeleteAttribute(string name)
        {
            if (H5A.exists(Id, name))
            {
                H5A.delete(Id, name);
                return true;
            }
            return false;
        }
        */
        #endregion

        protected override void WriteLow(HDFDatatype memtype, HDFDataspace memdspace, HDFDataspace filedspace, IntPtr buf)
        {
            H5D.write(Id, memtype.Id, (memdspace != null) ? memdspace.Id : H5S.ALL, (filedspace == null) ? H5S.ALL : filedspace.Id, 0, buf);
        }

        protected override void ReadLow(HDFDatatype memtype, HDFDataspace memdspace, HDFDataspace filedspace, IntPtr rawBytes)
        {
            H5D.read(Id, memtype.Id, (memdspace != null) ? memdspace.Id : H5S.ALL, (filedspace == null) ? H5S.ALL : filedspace.Id, 0, rawBytes);
        }

        #region IDataset Members


        public IReference CreateReference(IDataspace dataspace)
        {
            return new HDFReference(this, dataspace);
        }

        #endregion

        public long[] ChunkDimensions
        {
            get
            {
                var id = H5D.get_create_plist(Id);
                using (var plist = new HDFDatasetCreateProperty(id))
                using (var ds = Dataspace)
                    return plist.GetChunking(ds.Dimensions.Length);
            }
        }
    }


    public abstract class HDFGroupContainer : HDFIntraFile, IGroup
    {
        internal HDFGroupContainer(IChunkFile parent)
            : base(parent)
        {
        }

        public IDataset CreateDataset(string datasetName, IDatatype datatype, IDataspace dataspace, HDFDatasetCreateProperty prop)
        {
            return new HDFDataset(this, datasetName, datatype, dataspace, prop);
        }

        #region IContainer Members

        public IGroup CreateGroup(string groupName)
        {
            // If the group already exists, just return it
            IGroup g = (IGroup) GetChild(groupName);
            if(g != null)
                return g;
            else
                return new HDFGroup(this, groupName);
        }

        public IDataset CreateDataset(string datasetName, IDatatype datatype, IDataspace dataspace)
        {
            return new HDFDataset(this, datasetName, datatype, dataspace);
        }


        public IDataset OpenDataset(string datasetName, HDFDatasetAccessProperty plist)
        {
            return new HDFDataset(this, datasetName, plist);
        }

        public IChunkElement GetChild(string datasetName)
        {
            int childId = H5O.openNoThrow(Id, datasetName);

            if (childId < 0)
                return null;

            H5I.H5I_type type = H5I.get_type(childId);
            //H5O.H5O_type_t type = H5O.get_info_by_name(Id, datasetName).type;

            return GetByType(childId, type);
        }

        public bool RemoveChild(string datasetName)
        {
            if(H5L.exists(Id, datasetName))
            {
                H5L.delete(Id, datasetName);
                return true;
            }
            return false;
        }


        IChunkElement GetByType(int childId, H5O.H5O_type_t type)
        {
            switch (type)
            {
                case H5O.H5O_type_t.H5O_TYPE_DATASET:
                    return new HDFDataset(File, childId);

                case H5O.H5O_type_t.H5O_TYPE_GROUP:
                    return new HDFGroup(File, childId);
            }

            throw new ApplicationException("Unknown HDF child type" + type);  
        }


        IChunkElement GetByType(int childId, H5I.H5I_type type)
        {
            switch (type)
            {
                case H5I.H5I_type.H5I_DATASET:
                    return new HDFDataset(File, childId);

                case H5I.H5I_type.H5I_GROUP:
                    return new HDFGroup(File, childId);
            }

            throw new ApplicationException("Unknown HDF child type" + type);
        }


        IChunkElement GetChild(uint index)
        {
            H5O.H5O_type_t type = H5O.get_info_by_idx(Id, ".", index).type;

            int childId = H5O.open_by_idx(Id, ".", index);

            return GetByType(childId, type);
        }


        public IAttributeTarget GetNode(string datasetName)
        {
            return (IAttributeTarget)GetChild(datasetName);
        }


        public IChunkElement[] GetChildren()
        {
            ulong numchild = H5G.get_num_objs(Id);

            IChunkElement [] results = new IChunkElement[numchild];

            for (uint i = 0; i < numchild; i++)
                results[i] = GetChild(i);

            return results;
        }

        #endregion
    }

    /// <summary>
    /// Top level container for an HDF file
    /// </summary>
    public class HDFFile : HDFGroupContainer, IChunkFile
    {
        public string FileName { get; private set; }

        
        public static IChunkFile Open(string name, FileMode forWrite, FileAccess access)
        {
            return new HDFFile(name, forWrite, access);
        }


        public static IChunkFile Open(string name, FileMode forWrite, FileAccess access, HDFProperty accessPlist)
        {
            return new HDFFile(name, forWrite, access, accessPlist);
        }



        private HDFFile(string filename, FileMode forWrite, FileAccess access) : this(filename, forWrite, access, null)
        {

        }

        private HDFFile(string filename, FileMode forWrite, FileAccess access, HDFProperty fileAccessProperties) : base(null)
        {
            FileName = filename;

            // We want error messages
            H5E.ShowErrors = true;

            H5F.AccessMode h5Access = (access == FileAccess.ReadWrite) ? H5F.AccessMode.ACC_RDWR : H5F.AccessMode.ACC_RDONLY;
            
            if(forWrite == FileMode.OpenOrCreate)
                forWrite = System.IO.File.Exists(filename) ? FileMode.Open : FileMode.Create;

            switch(forWrite)
            {
                case FileMode.Open:
                    // Check ourselves -- will give a better error message.
                    if (!System.IO.File.Exists(filename))
                        throw new FileNotFoundException(filename);
                    var fileAccessId = fileAccessProperties != null ? fileAccessProperties.Id : H5P.DEFAULT; 
                    Id = H5F.open(filename, h5Access, fileAccessId);
                    break;
                case FileMode.Create:
                    Id = H5F.create(filename, H5F.AccessMode.ACC_TRUNC);
                    break;
                default:
                    throw new NotSupportedException("Unsupported access mode for HDF files");

            }
        }


        /// <summary>
        /// files can't have attributes - just return an empty array
        /// </summary>
        /// <returns></returns>
        public override IDataContainer[] GetAttributes()
        {
            return new IDataContainer[0];
        }
        
        #region IChunkElement Members

        public override string Name
        {
            get { return FileName; }
        }

        #endregion

        #region IDisposable Members

        protected override void Dispose(bool disposing)
        {
            if (Id != -1)
            {
                H5F.close(Id);
                Id = -1;
            }
        }

        #endregion

        #region IChunkFile Members

        public IDatatype CreateDatatype(Type t)
        {
            // If the user passes in an array type we assume that the dataspace we want is actually typed by the ELEMENTS of the array
            if (t.IsArray)
                t = t.GetElementType();

            return new HDFDatatype(t);
        }

        #endregion


        #region IChunkFile Members


        public IDataspace CreateDataspace(long[] curDimensions, long[] maxDimensions)
        {
            return new HDFDataspace(File, curDimensions, maxDimensions);
        }

        public IDataspace CreateDataspace()
        {
            return new HDFDataspace(File);
        }



        public override IChunkFile File
        {
            get
            {
                return this;
            }
        }

        #endregion


        #region IChunkFile Members

        private bool flatten = false;

        public bool FlattenArrays
        {
            get
            {
                return flatten;
            }
            set
            {
                flatten = value;
            }
        }

        #endregion
    }
}
