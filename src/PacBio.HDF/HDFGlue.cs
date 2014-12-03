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
using System.IO;
using System.Text;
using System.Runtime.InteropServices;
using System.Diagnostics;

// Don't worry about the HDF function names
// ReSharper disable InconsistentNaming

namespace PacBio.HDF
{
    [Serializable]
    public class NotFoundException : Exception
    {
    }


    public class HDFGlue
    {
        // Calls to the HDF5 library should be serialized through this lock.
        // This prevents us from having to build the pthread thread-safe version of 
        // which apparently no longer is available on windows.  The thread-safe HDF5 just serializes the 
        // calls anyway, so maybe the runtime can do something smarter if we handle the locking
        // at the VM level?
        public static readonly object HDFLock = new Object();
        
        // This is the string that get used in P/Invoke
        internal const string hdfPInvokeName = "hdf5";
        internal const string hdfHlPInvokeName = "hdf5_hl";

        internal const string hdfUnix = "libhdf5.so";
        internal const string hdfWindows = "hdf5.dll";
        internal const string hdfMac = "libhdf5.dylib";


        private static string HDFLibName
        {
            get
            {
                if (CRuntime.RunningPlatform == CRuntime.Platform.Windows)
                {
                    return hdfWindows;
                }
                else if (CRuntime.RunningPlatform == CRuntime.Platform.Mac)
                {
                    return hdfMac;
                }
                else
                {
                    return hdfUnix;
                }
            }
        }

        private static IntPtr dllptr = IntPtr.Zero;

        /// <summary>
        /// HDF hides nasty ints in _g globals - this is how we find em
        /// </summary>
        /// <param name="globalname"></param>
        /// <returns></returns>
        internal static int GetHDFGlobal(string globalname)
        {
            if (dllptr.Equals(IntPtr.Zero))
            {
                // Try for the hdf5 dll packaged with the assembly
                dllptr = CRuntime.GetModuleHandle(HDFLibName);

                if (dllptr == IntPtr.Zero)
                {
                    throw new ApplicationException(String.Format("Can't find HDF library: {0}", HDFLibName));
                }
            }

            IntPtr addr = CRuntime.GetProcAddress(dllptr, globalname);

            if (addr == IntPtr.Zero)
            {
                throw new Exception("Couldn't find global " + globalname + " in dll " + HDFLibName);
            }

            int id = Marshal.ReadInt32(addr);
            return id;
        }

        /// <summary>
        /// A flag for checking that HDF5 is initialized properly.  Set to true after a successful call to H5open()
        /// In the constructor of HDFGlue.
        /// </summary>
        public static bool Available { get; private set; }

        /// <summary>
        /// Static constructor to ensure HDF library is loaded
        /// </summary>
        static HDFGlue()
        {
            Available = false;

            lock (HDFLock)
            {
                CheckError(H5open());
            }

            Available = true;
        }

        protected delegate int HDFCall();

        /// <summary>
        /// Check for an HDF error and assume the user is trying to open an item
        /// </summary>
        /// May throw NotFoundException
        protected static int CheckOpenError(HDFCall funct)
        {
            bool oldErr = H5E.ShowErrors;
            H5E.ShowErrors = false;

            int errcode = funct();

            H5E.ShowErrors = oldErr;

            if (errcode < 0)
            {
                throw new NotFoundException();
            }

            return errcode;
        }


        /// <summary>
        /// Check for an HDF error and assume the user is trying to open an item
        /// </summary>
        /// May throw NotFoundException
        protected static int CheckOpenNoThrow(HDFCall funct)
        {
            bool oldErr = H5E.ShowErrors;
            H5E.ShowErrors = false;

            int errcode = funct();

            H5E.ShowErrors = oldErr;

            return errcode;
        }


        /// <summary>
        /// Check for an HDF error and throw an exception if a problem is encountered
        /// </summary>
        /// <param name="errcode"></param>
        /// <returns></returns>
        protected static int CheckError(int errcode)
        {
            if (errcode < 0)
            {
                throw new IOException("HDF5 error=" + errcode);
            }

            return errcode;
        }

        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern int H5open();
    }

    /// <summary>
    /// pinvokes to call the c style HDF library
    /// </summary>
    public class H5A : HDFGlue
    {
        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern int H5Acreate2(int loc, string name, int type_id, int space_id, int acpl_id, int aapl_id);

        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern int H5Aopen(int id,
                                          string name,
                                          int accessplist);

        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern int H5Aexists(int id,
                                            string name);

        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern int H5Aclose(int id);

        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern int H5Adelete(int id, string name);

        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern int H5Aget_type(int id);

        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern int H5Aget_space(int id);

        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern int H5Aget_num_attrs(int id);

        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern int H5Aopen_idx(int loc_id, uint idx);

        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern int H5Awrite(int id, int memtypeid, IntPtr buf);

        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern int H5Aread(int dataset_id, int mem_type_id, [Out] IntPtr buf);

        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern int H5Aget_name(int attr_id, int buf_size, StringBuilder buf);

        public static int create(int fileId, string name, int datatype, int dataspace)
        {
            lock (HDFLock)
            {
                return CheckError(H5Acreate2(fileId, name, datatype, dataspace, H5P.DEFAULT, H5P.DEFAULT));
            }
        }

        public static int open(int fileId, string name)
        {
            lock (HDFLock)
            {
                return CheckError(H5Aopen(fileId, name, H5P.DEFAULT));
            }
        }

        public static int delete(int fileId, string name)
        {
            lock (HDFLock)
            {
                return CheckError(H5Adelete(fileId, name));
            }
        }

        /// <summary>
        /// true if the named property exists on the specified object
        /// </summary>
        /// <param name="fileId"></param>
        /// <param name="name"></param>
        /// <returns></returns>
        public static bool exists(int fileId, string name)
        {
            lock (HDFLock)
            {
                return CheckError(H5Aexists(fileId, name)) > 0;
            }
        }

        public static int write(int id, int memtypeid, IntPtr obj)
        {
            lock (HDFLock)
            {
                return CheckError(H5Awrite(id, memtypeid, obj));
            }
        }

        public static int read(int id, int memtypeid, IntPtr obj)
        {
            lock (HDFLock)
            {
                return CheckError(H5Aread(id, memtypeid, obj));
            }
        }

        public static string get_name(int id)
        {
            StringBuilder builder;

            lock (HDFLock)
            {
                int sz = CheckError(H5Aget_name(id, 0, null));

                builder = new StringBuilder(sz + 1);
                CheckError(H5Aget_name(id, builder.Capacity, builder));
            }

            return builder.ToString();
        }

        public static int get_num_attrs(int id)
        {
            lock (HDFLock)
            {
                return CheckError(H5Aget_num_attrs(id));
            }
        }

        public static int open_idx(int id, uint index)
        {
            lock (HDFLock)
            {
                return CheckError(H5Aopen_idx(id, index));
            }
        }

        public static int get_type(int id)
        {
            lock (HDFLock)
            {
                return CheckError(H5Aget_type(id));
            }
        }

        public static int get_space(int id)
        {
            lock (HDFLock)
            {
                return CheckError(H5Aget_space(id));
            }
        }

        public static void close(int id)
        {
            lock (HDFLock)
            {
                CheckError(H5Aclose(id));
            }
        }
    }

    public class H5F : HDFGlue
    {
        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern int H5Fcreate(string name,
                                            uint flags,
                                            int createplist,
                                            int accessplist);

        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern int H5Fopen(string name,
                                          uint flags,
                                          int accessplist);

        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern int H5Fclose(int id);

        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern int H5Fget_obj_count(int file_id, uint types);

        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern int H5Fget_obj_ids(int file_id, uint types, int maxids, [Out] int[] ids);


        [Flags]
        public enum AccessMode
        {
            ACC_RDONLY = 0x0000, /*absence of rdwr => rd-only */
            ACC_RDWR = 0x0001, /*open for read and write    */
            ACC_TRUNC = 0x0002, /*overwrite existing files   */
            ACC_EXCL = 0x0004, /*fail if file already exists*/
            ACC_DEBUG = 0x0008, /*print debug info	        */
            ACC_CREAT = 0x0010 /*create non-existing files  */
        }

        [Flags]
        public enum ObjTypes
        {
            OBJ_FILE = (0x0001), /* File objects */
            OBJ_DATASET = (0x0002), /* Dataset objects */
            OBJ_GROUP = (0x0004), /* Group objects */
            OBJ_DATATYPE = (0x0008), /* Named datatype objects */
            OBJ_ATTR = (0x0010), /* Attribute objects */
            OBJ_ALL = (OBJ_FILE | OBJ_DATASET | OBJ_GROUP | OBJ_DATATYPE | OBJ_ATTR),
            OBJ_LOCAL = (0x0020) /* Restrict search to objects opened through current file ID */
        }

        public static int create(string name, AccessMode createMode)
        {
            lock (HDFLock)
            {
                return CheckError(H5Fcreate(name, (uint) createMode, H5P.DEFAULT, H5P.DEFAULT));
            }
        }

        public static int open(string name, AccessMode mode)
        {
            lock (HDFLock)
            {
                return CheckError(H5Fopen(name, (uint) mode, H5P.DEFAULT));
            }
        }

        public static int open(string name, AccessMode mode, int accessPlist)
        {
            lock (HDFLock)
            {
                return CheckError(H5Fopen(name, (uint) mode, accessPlist));
            }
        }

        public static int close(int id)
        {
            lock (HDFLock)
            {
                return CheckError(H5Fclose(id));
            }
        }

        public static int get_obj_count(int id, ObjTypes types)
        {
            lock (HDFLock)
            {
                return CheckError(H5Fget_obj_count(id, (uint) types));
            }
        }

        public static int[] get_obj_ids(int id, ObjTypes types)
        {
            var ids = new int[get_obj_count(id, types)];

            lock (HDFLock)
            {
                if (ids.Length != 0)
                    CheckError(H5Fget_obj_ids(id, (uint) types, ids.Length, ids));
            }

            return ids;
        }
    }

    public class H5D : HDFGlue
    {
        /// <summary>
        /// Chunking optiosn for use with H5Pset_layout
        /// </summary>
        public int COMPACT = 0; /*raw data is very small		     */

        public int CONTIGUOUS = 1; /*the default				     */
        public int CHUNKED = 2; /*slow and fancy			     */

        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern int H5Dcreate2(int fileid,
                                             string name,
                                             int dataType,
                                             int dataSpace,
                                             int linkprops,
                                             int createprops,
                                             int accessprops);


        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern int H5Dopen2(int id,
                                           string name,
                                           int accessplist);


        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern int H5Dclose(int id);

        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern int H5Dget_type(int id);

        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern int H5Dget_space(int id);

        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern int H5Dget_create_plist(int id);


        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern int H5Dextend(int space_id, [In] long[] dims);

        /// <summary>
        /// HDF data write
        /// </summary>
        /// <param name="id"></param>
        /// <param name="memtypeid"></param>
        /// <param name="memspaceid"></param>
        /// <param name="filespaceid"></param>
        /// <param name="xferplist"></param>
        /// <param name="buf"></param>
        /// <returns></returns>
        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern int H5Dwrite(int id, int memtypeid, int memspaceid, int filespaceid, int xferplist,
                                           IntPtr buf);

        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern int H5Dread(int dataset_id, int mem_type_id, int mem_space_id, int file_space_id,
                                          int xfer_plist_id, [Out] IntPtr buf);

        public static int create(int fileId, string name, int datatype, int dataspace, int linkPropId,
                                 int creationPropId)
        {
            lock (HDFLock)
            {
                return CheckError(H5Dcreate2(fileId, name, datatype, dataspace, linkPropId, creationPropId, H5P.DEFAULT));
            }
        }

        public static int open(int fileId, string name)
        {
            lock (HDFLock)
            {
                return CheckError(H5Dopen2(fileId, name, H5P.DEFAULT));
            }
        }

        public static int open(int fileId, string name, int accessPlist)
        {
            lock (HDFLock)
            {
                return CheckError(H5Dopen2(fileId, name, accessPlist));
            }
        }

        public static int write(int id, int memtypeid, int memspaceid, int filespaceid, int xferplist, IntPtr obj)
        {
            lock (HDFLock)
            {
                return CheckError(H5Dwrite(id, memtypeid, memspaceid, filespaceid, xferplist, obj));
            }
        }

        public static int read(int id, int memtypeid, int memspaceid, int filespaceid, int xferplist, IntPtr obj)
        {
            lock (HDFLock)
            {
                return CheckError(H5Dread(id, memtypeid, memspaceid, filespaceid, xferplist, obj));
            }
        }


        public static int extend(int id, long[] dims)
        {
            lock (HDFLock)
            {
                return CheckError(H5Dextend(id, dims));
            }
        }


        public static int get_type(int id)
        {
            lock (HDFLock)
            {
                return CheckError(H5Dget_type(id));
            }
        }

        public static int get_space(int id)
        {
            lock (HDFLock)
            {
                return CheckError(H5Dget_space(id));
            }
        }

        public static int get_create_plist(int id)
        {
            lock (HDFLock)
            {
                return CheckError(H5Dget_create_plist(id));
            }
        }

        public static void close(int id)
        {
            lock (HDFLock)
            {
                CheckError(H5Dclose(id));
            }
        }
    }


    public class H5O : HDFGlue
    {
        private const int H5_INDEX_NAME = 0; /* Index on names 			*/
        private const int H5_INDEX_CRT_ORDER = 0; /* Index on creation order 		*/

        private const int H5_ITER_INC = 0,
                          /* Increasing order */
                          H5_ITER_DEC = 1,
                          /* Decreasing order */
                          H5_ITER_NATIVE = 2; /* No particular order, whatever is fastest */

        public enum H5O_type_t
        {
            H5O_TYPE_UNKNOWN = -1, /* Unknown object type		*/
            H5O_TYPE_GROUP, /* Object is a group		*/
            H5O_TYPE_DATASET, /* Object is a dataset		*/
            H5O_TYPE_NAMED_DATATYPE, /* Object is a named data type	*/
            H5O_TYPE_NTYPES /* Number of different object types (must be last!) */
        }

        [StructLayout(LayoutKind.Sequential)]
        public struct H5_ih_info_t
        {
            public ulong index_size; /* btree and/or list */
            public ulong heap_size;
        }

        /* Information struct for object (for H5Oget_info/H5Oget_info_by_name/H5Oget_info_by_idx) */

        [StructLayout(LayoutKind.Sequential)]
        public struct H5O_info_t
        {
            public uint fileno; /* (unsigned long) File number that object is located in */
            public ulong addr; /* Object address in file	*/
            public H5O_type_t type; /* Basic object type (group, dataset, etc.) */
            public uint rc; /* Reference count of object    */

            public ulong atime; /* Access time			*/
            public ulong mtime; /* Modification time		*/
            public ulong ctime; /* Change time			*/
            public ulong btime; /* Birth time			*/
            public ulong num_attrs; /* # of attributes attached to object */

            // struct {
            public uint version; /* Version number of header format in file */
            public uint nmesgs; /* Number of object header messages */
            public uint nchunks; /* Number of object header chunks */
            public uint flags; /* Object header status flags */
            // struct {
            public ulong total; /* Total space for storing object header in file */
            public ulong meta; /* Space within header for object header metadata information */
            public ulong mesg; /* Space within header for actual message information */
            public ulong free; /* Free space within object header */
            // } space;
            // struct {
            public ulong present; /* Flags to indicate presence of message type in header */
            public ulong shared; /* Flags to indicate message type is shared in header */
            // } mesg;
            // } hdr;
            /* Extra metadata storage for obj & attributes */
            // struct {
            public H5_ih_info_t obj; /* v1/v2 B-tree & local/fractal heap for groups, B-tree for chunked datasets */
            public H5_ih_info_t attr; /* v2 B-tree & heap for attributes */
            // } meta_size;
        }

        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern int H5Oopen_by_idx(int loc_id, string group_name, int index_field, int order, ulong n,
                                                 int lapl_id);

        // H5Oget_info_by_idx(hid_t loc_id, const char *group_name,
        // H5_index_t idx_type, H5_iter_order_t order, hsize_t n, H5O_info_t *oinfo, hid_t lapl_id);
        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern int H5Oget_info_by_idx(int loc_id, string group_name, int index_type, int order,
                                                     ulong n, out H5O_info_t object_info, int lapl_id);

        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern int H5Oget_info_by_name(int loc_id, string object_name,
                                                      out H5O_info_t object_info, int lapl_id);

        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern int H5Oopen(int loc_id, string object_name, int lapl_id);

        public static int open_by_idx(int loc_id, string group_name, ulong n)
        {
            lock (HDFLock)
            {
                return CheckError(H5Oopen_by_idx(loc_id, group_name, H5_INDEX_CRT_ORDER, H5_ITER_INC, n, H5P.DEFAULT));
            }
        }

        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern int H5Ocopy(int src_loc_id, string src_name, int dst_loc_id, string dst_name,
                                          int ocpypl_id, int lcpl_id);

        public static int open(int loc_id, string group_name)
        {
            lock (HDFLock)
            {
                return CheckOpenError(() => H5Oopen(loc_id, group_name, H5P.DEFAULT));
            }
        }

        public static int openNoThrow(int loc_id, string group_name)
        {
            lock (HDFLock)
            {
                return CheckOpenNoThrow(() => H5Oopen(loc_id, group_name, H5P.DEFAULT));
            }
        }

        public static H5O_info_t get_info_by_idx(int loc_id, string group_name, ulong n)
        {
            H5O_info_t result = new H5O_info_t();
            lock (HDFLock)
            {
                CheckError(H5Oget_info_by_idx(loc_id, group_name, H5_INDEX_CRT_ORDER, H5_ITER_INC, n, out result,
                                              H5P.DEFAULT));
            }
            return result;
        }

        public static H5O_info_t get_info_by_name(int loc_id, string name)
        {
            H5O_info_t result = new H5O_info_t();

            lock (HDFLock)
            {
                CheckError(H5Oget_info_by_name(loc_id, name, out result, H5P.DEFAULT));
            }

            return result;
        }

        public static int copy(int src_loc_id, string src_name, int dest_loc_id, string dest_name)
        {
            lock (HDFLock)
            {
                return CheckError(H5Ocopy(src_loc_id, src_name, dest_loc_id, dest_name, 0, 0));
            }
        }
    }


    public class H5G : HDFGlue
    {
        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern int H5Gcreate2(int fileId,
                                             string groupName,
                                             int linkCreation,
                                             int groupCreation,
                                             int groupAccess);

        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern int H5Gclose(int fileId);

        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern int H5Gget_objtype_by_idx(int file_id, uint idx);

        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern int H5Gget_num_objs(int loc_id, out ulong num_obj);

        public static int create(int fileId, string groupName)
        {
            lock (HDFLock)
            {
                return CheckError(H5Gcreate2(fileId, groupName, H5P.DEFAULT, H5P.DEFAULT, H5P.DEFAULT));
            }
        }

        public static void close(int id)
        {
            lock (HDFLock)
            {
                CheckError(H5Gclose(id));
            }
        }

        public static ulong get_num_objs(int id)
        {
            ulong num;

            lock (HDFLock)
            {
                CheckError(H5Gget_num_objs(id, out num));
            }

            return num;
        }
    }

    /// <summary>
    /// Find info about HDF5 identifiers
    /// </summary>
    public class H5I : HDFGlue
    {
        public enum H5I_type
        {
            H5I_FILE = 1, /*type ID for File objects		    */
            H5I_GROUP, /*type ID for Group objects		    */
            H5I_DATATYPE, /*type ID for Datatype objects		    */
            H5I_DATASPACE, /*type ID for Dataspace objects		    */
            H5I_DATASET, /*type ID for Dataset objects		    */
            H5I_ATTR, /*type ID for Attribute objects		    */
            H5I_REFERENCE, /*type ID for Reference objects		    */
            H5I_VFL, /*type ID for virtual file layer	    */
            H5I_GENPROP_CLS, /*type ID for generic property list classes */
            H5I_GENPROP_LST, /*type ID for generic property lists        */
            H5I_ERROR_CLASS, /*type ID for error classes		    */
            H5I_ERROR_MSG, /*type ID for error messages		    */
            H5I_ERROR_STACK, /*type ID for error stacks		    */
        }

        /// <summary>
        /// Used to iterate over object tables
        /// </summary>
        /// <param name="objptr"></param>
        /// <param name="id"></param>
        /// <param name="userData"></param>
        /// <returns>zero to keep searching</returns>
        public delegate int SearchCallback(IntPtr objptr, int id, IntPtr userData);

        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern IntPtr H5Isearch(int type, SearchCallback func, IntPtr userData);

        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern int H5Iget_name(int obj_id, StringBuilder name, int size);

        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern int H5Iget_type(int obj_id);

        public static string get_name(int objid)
        {
            StringBuilder build = new StringBuilder(256);

            lock (HDFLock)
            {
                CheckError(H5Iget_name(objid, build, build.Capacity));
            }

            return build.ToString();
        }

        public static H5I_type get_type(int obj_id)
        {
            int code;
            lock (HDFLock)
            {
                code = CheckError(H5Iget_type(obj_id));
            }

            return (H5I_type) code;
        }


        public static IntPtr search(H5I_type type, SearchCallback callback, IntPtr userData)
        {
            lock (HDFLock)
            {
                return H5Isearch((int) type, callback, userData);
            }
        }
    }

    public class H5P : HDFGlue
    {
        /// <summary>
        /// H5P_DEFAULT
        /// </summary>
        public const int DEFAULT = 0;

        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern int H5Pset_layout(int plist, int layout);

        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern int H5Pcreate(int plist_class);

        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern int H5Pclose(int pid);

        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern int H5Pset_chunk(int plist, int ndims, [In] long[] dim);

        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern int H5Pget_chunk(int plist, int ndims, [Out] long[] dim);


        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern int H5Pset_cache(int plist, int mdc_nelmts, int rdcc_nelmts, ulong rdcc_bytes,
                                               double rdcc_w0);

        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern int H5Pset_chunk_cache(int plist, IntPtr rdcc_nslots, IntPtr rdcc_nbytes, double rdcc_w0);

        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern int H5Pset_deflate(int plits, uint level);


        //herr_t H5Pset_chunk_cache( hid_t dapl_id, size_t rdcc_nslots, size_t rdcc_nbytes, double rdcc_w0 )

        /// <summary>
        /// Get the name of this property
        /// </summary>
        /// <param name="pcid"></param>
        /// <returns></returns>
        /// This method _will_ cause a memory leak, but should only be calling it a few times at startup
        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern string H5Pget_class_name(int pcid);


        public static int create(string prop_class_id)
        {
            // Need to call this to ensure that the HDF5 dll has been loaded
            // into the process, otherwise GetModuleHandle below will return null
            H5E.ShowErrors = true;

            // We find our int by looking in the hdf5 DLL
            int propid = HDFGlue.GetHDFGlobal(prop_class_id + "_g");

            lock (HDFLock)
            {
                return CheckError(H5Pcreate(propid));
            }
        }

        public static void set_chunk(int id, long[] dims)
        {
            lock (HDFLock)
            {
                CheckError(H5Pset_chunk(id, dims.Length, dims));
            }
        }

        public static void set_deflate(int id, uint level)
        {
            lock (HDFLock)
            {
                CheckError(H5Pset_deflate(id, level));
            }
        }

        public static long[] get_chunk(int id, int ndims)
        {
            var chunking = new long[ndims];
            lock (HDFLock)
            {
                CheckError(H5Pget_chunk(id, ndims, chunking));
            }

            return chunking;
        }

        public static void set_cache(int plistId, int mdc_nelmts, int rdcc_nelmts, ulong rdcc_bytes,
                                     double rdcc_w0)
        {
            lock (HDFLock)
            {
                CheckError(H5Pset_cache(plistId, mdc_nelmts, rdcc_nelmts, rdcc_bytes, rdcc_w0));
            }
        }

        public static void set_chunk_cache(int plistId, uint rdcc_nslots, uint rdcc_nbytes, double rdcc_w0)
        {
            lock (HDFLock)
            {
                CheckError(H5Pset_chunk_cache(plistId, (IntPtr) rdcc_nslots, (IntPtr) rdcc_nbytes, rdcc_w0));
            }
        }


        public static int set_layout(int id, int layout)
        {
            lock (HDFLock)
            {
                return CheckError(H5Pset_layout(id, layout));
            }
        }

        public static int close(int id)
        {
            lock (HDFLock)
            {
                return CheckError(H5Pclose(id));
            }
        }
    }


    /// <summary>
    /// HDF library exceptions wrapper
    /// </summary>
    public class HDFException : IOException
    {
        public HDFException(string msg)
            : base(msg)
        {
        }
        public HDFException(string msg, Exception inner) :base(msg, inner){}

        public string Description;
        public string MajorMessage;
        public string MinorMessage;
    }


    /// <summary>
    /// HDF Error printing/parsing
    /// </summary>
    public class H5E : HDFGlue
    {
        // typedef herr_t (*H5E_walk_t)(int n, H5E_error_t *err_desc, void *client_data) 

        /*
         * typedef struct H5E_error2_t {
        hid_t       cls_id;         
        hid_t       maj_num;	
        hid_t       min_num;	
        unsigned	line;		
        const char	*func_name;   	
        const char	*file_name;	
        const char	*desc;		
        } H5E_error2_t;
        */

        [StructLayout(LayoutKind.Sequential)]
        private struct error_t
        {
            public int cls_id,
                       maj_num,
                       min_num;

            public uint line;

            [MarshalAs(UnmanagedType.LPStr)] public string func_name,
                                                           file_name,
                                                           desc;
        }


        //MarshalAs(UnmanagedType.LPStruct)
        private delegate int ErrorCallback(int n, [In] ref error_t errDesc, IntPtr userData);

        /// <summary>
        /// A hook for printing HDF errors
        /// </summary>
        /// <param name="whichstack"></param>
        /// <param name="cstream">A c style FILE * (useless to us)</param>
        /// <returns></returns>
        private delegate int PrintCallback(int whichstack, IntPtr cstream);

        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern int H5Ewalk2(int whichstack, int direction, ErrorCallback func, IntPtr userData);

        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern int H5Eset_auto2(int whichstack, PrintCallback func, IntPtr userData);

        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern int H5Eget_msg(int mesg_id, out int mesg_type, StringBuilder mesg, uint size);

        private const int WALK_UPWARD = 0; /*begin deep, end at API function    */
        private const int WALK_DOWNWARD = 1; /*begin at API function, end deep    */

        private static int MyErrorCallback(int n, [In] ref error_t errDesc, IntPtr userData)
        {
            var majorMessage = get_msg(errDesc.maj_num);
            var minorMessage = get_msg(errDesc.min_num);
            var description = errDesc.desc;

            string str = n + ": " + errDesc.file_name + ":" + errDesc.func_name + ":" +
                         errDesc.line + " " + majorMessage + "/" + minorMessage + ": " +
                         description;

            throw new HDFException("HDF: " + str)
                      {

                          MajorMessage = majorMessage,
                          MinorMessage = minorMessage,
                          Description = description
                      };

            // We just throw now, because something is hosed in the C HDF stack by this point
            // return 0;
        }

        private static int MyPrintCallback(int whichstack, IntPtr cstream)
        {
            return walk(whichstack);
        }

        private static string get_msg(int msgId)
        {
            StringBuilder build = new StringBuilder(256);
            int typ;

            lock (HDFLock)
            {
                CheckError(H5Eget_msg(msgId, out typ, build, (uint) build.Capacity));
            }

            return build.ToString();
        }


        /// <summary>
        /// Walk the error stack and print it
        /// </summary>
        /// <returns></returns>
        private static int walk(int whichstack)
        {
            ErrorCallback cb = new ErrorCallback(MyErrorCallback);
            int result;

            // only stack 0 is defined right now
            lock (HDFLock)
            {
                result = H5Ewalk2(whichstack, WALK_UPWARD, cb, IntPtr.Zero);
            }

            return result;
        }

        private static PrintCallback cb = new PrintCallback(MyPrintCallback);

        private static bool showErr = false;

        public static bool ShowErrors
        {
            get { return showErr; }

            set
            {
                showErr = value;
                lock (HDFLock)
                {
                    if (showErr)
                        CheckError(H5Eset_auto2(0, cb, IntPtr.Zero));
                    else
                        CheckError(H5Eset_auto2(0, null, IntPtr.Zero));
                }
            }
        }
    }


    /// <summary>
    /// HDF5 dataspace API
    /// </summary>
    public class H5S : HDFGlue
    {
        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern int H5Screate(int dataspacetype);

        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern int H5Screate_simple(int rank, [In] long[] dims, [In] long[] maxdims);

        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern int H5Sclose(int id);

        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern int H5Scopy(int id);

        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern int H5Sset_extent_simple(int space_id, int rank, [In] long[] curdims, [In] long[] maxdims);

        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern int H5Sget_simple_extent_ndims(int space_id);

        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern long H5Sget_select_npoints(int space_id);

        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern int H5Sget_simple_extent_dims(int space_id, [Out] long[] dims, [Out] long[] maxdims);

        private enum seloper_t
        {
            H5S_SELECT_NOOP = -1, /* error                                     */
            H5S_SELECT_SET = 0, /* Select "set" operation 		     */
            H5S_SELECT_OR,
            /* Binary "or" operation for hyperslabs
                                         * (add new selection to existing selection)
                                         * Original region:  AAAAAAAAAA
                                         * New region:             BBBBBBBBBB
                                         * A or B:           CCCCCCCCCCCCCCCC
                                         */
            H5S_SELECT_AND,
            /* Binary "and" operation for hyperslabs
                                         * (only leave overlapped regions in selection)
                                         * Original region:  AAAAAAAAAA
                                         * New region:             BBBBBBBBBB
                                         * A and B:                CCCC
                                         */
            H5S_SELECT_XOR,
            /* Binary "xor" operation for hyperslabs
                                         * (only leave non-overlapped regions in selection)
                                         * Original region:  AAAAAAAAAA
                                         * New region:             BBBBBBBBBB
                                         * A xor B:          CCCCCC    CCCCCC
                                         */
            H5S_SELECT_NOTB,
            /* Binary "not" operation for hyperslabs
                                         * (only leave non-overlapped regions in original selection)
                                         * Original region:  AAAAAAAAAA
                                         * New region:             BBBBBBBBBB
                                         * A not B:          CCCCCC
                                         */
            H5S_SELECT_NOTA,
            /* Binary "not" operation for hyperslabs
                                         * (only leave non-overlapped regions in new selection)
                                         * Original region:  AAAAAAAAAA
                                         * New region:             BBBBBBBBBB
                                         * B not A:                    CCCCCC
                                         */
            H5S_SELECT_APPEND, /* Append elements to end of point selection */
            H5S_SELECT_PREPEND, /* Prepend elements to beginning of point selection */
            H5S_SELECT_INVALID /* Invalid upper bound on selection operations */
        }

        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern int H5Sselect_hyperslab(int space_id, seloper_t op,
                                                      [In] long[] start,
                                                      [In] long[] stride,
                                                      [In] long[] count,
                                                      [In] long[] block);

        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern int H5Sselect_all(int space_id);

        /// <summary>
        /// Dataspace types
        /// </summary>
        public const int SCALAR = 0; /*scalar variable                            */

        public const int SIMPLE = 1; /*simple data space                          */
        public const int NULL = 2; /*null data space                            */


        /// <summary>
        /// A simple in memory H5S_ALL dataset
        /// </summary>
        public const int ALL = 0;

        public static int create(int dspacetype)
        {
            lock (HDFLock)
            {
                return CheckError(H5Screate(dspacetype));
            }
        }

        public static int create_simple(long[] dims, long[] maxdims)
        {
            lock (HDFLock)
            {
                int rank = maxdims.Length;
                return CheckError(H5Screate_simple(rank, dims, maxdims));
            }
        }


        public static int copy(int dspacetype)
        {
            lock (HDFLock)
            {
                return CheckError(H5Scopy(dspacetype));
            }
        }


        public static long get_select_npoints(int dspacetype)
        {
            lock (HDFLock)
            {
                long res;
                CheckError((int) (res = H5Sget_select_npoints(dspacetype)));
                return res;
            }
        }


        public static void select_hyperslab(int id,
                                            long[] start,
                                            long[] stride,
                                            long[] count,
                                            long[] block)
        {
            lock (HDFLock)
            {
                CheckError(H5Sselect_hyperslab(id, seloper_t.H5S_SELECT_SET, start, stride, count, block));
            }
        }


        public static void select_all(int id)
        {
            lock (HDFLock)
            {
                CheckError(H5Sselect_all(id));
            }
        }


        public static void get_simple_extent(int id, out long[] curdims, out long[] maxdims)
        {
            lock (HDFLock)
            {
                int rank = CheckError(H5Sget_simple_extent_ndims(id));

                curdims = new long[rank];
                maxdims = new long[rank];

                // scalar dataspaces will have rank zero
                if (rank != 0)
                    CheckError(H5Sget_simple_extent_dims(id, curdims, maxdims));
            }
        }

        public static int set_extent_simple(int id, long[] curdims, long[] maxdims)
        {
            Debug.Assert(curdims.Length == maxdims.Length);

            lock (HDFLock)
            {
                return CheckError(H5Sset_extent_simple(id, curdims.Length, curdims, maxdims));
            }
        }

        public static void close(int id)
        {
            lock (HDFLock)
            {
                CheckError(H5Sclose(id));
            }
        }
    }

    /// <summary>
    /// The 'high level' HDF5 utility library
    /// </summary>
    public class H5L : HDFGlue
    {
        /// <summary>
        /// 
        /// </summary>
        /// <param name="datatypename"></param>
        /// <param name="langtype">Only langtype 0 is supprted at present</param>
        /// <returns></returns>
        [DllImport(hdfHlPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern int H5LTtext_to_dtype(string datatypename, int langtype);


        // Use a UIntPtr because on x86_64 Linux the final parameter is a pointer
        // to a 64-bit integer (it is a size_t*).
        [DllImport(hdfHlPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern int H5LTdtype_to_text(int datatype, StringBuilder outstr, int lang_type,
                                                    [In, Out] ref UIntPtr outstrlen);

        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern int H5Ldelete(int loc_id, string name, int lapl_id);

        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern int H5Lexists(int loc_id, string name, int lapl_id);


        /// <summary>
        /// Given a string containing a datatype name, return the integer data type ID
        /// </summary>
        /// <param name="datatypename"></param>
        /// <returns></returns>
        public static int text_to_dtype(string datatypename)
        {
            lock (HDFLock)
                return CheckError(H5LTtext_to_dtype(datatypename, 0));
        }

        public static string dtype_to_text(int dtype)
        {
            StringBuilder builder;
            // First find out how many bytes we will need
            lock (HDFLock)
            {
                UIntPtr outstrlen = UIntPtr.Zero;

                CheckError(H5LTdtype_to_text(dtype, null, 0, ref outstrlen));
                builder = new StringBuilder((int) outstrlen);
                CheckError(H5LTdtype_to_text(dtype, builder, 0, ref outstrlen));
            }
            return builder.ToString();
        }

        public static void delete(int loc_id, string name)
        {
            lock (HDFLock)
                CheckError(H5Ldelete(loc_id, name, H5P.DEFAULT));
        }

        public static bool exists(int loc_id, string name)
        {
            int r;

            lock (HDFLock)
                r = CheckError(H5Lexists(loc_id, name, H5P.DEFAULT));

            return r > 0 ? true : false;
        }
    }

    /// <summary>
    /// HDF5 reference support
    /// </summary>
    public class H5R : HDFGlue
    {
        [StructLayout(LayoutKind.Sequential, Pack = 4)]
        public struct RegRef
        {
            // 12 is a 'magic' translation from H5Rpublic.h/hdset_reg_ref_t
            [MarshalAs(UnmanagedType.ByValArray, ArraySubType = UnmanagedType.Struct, SizeConst = 12)] public byte[] ptr;
                                                                                                                     
        }

        // ref_types for H5Rcreate
        private const int H5R_OBJECT = 0;
        private const int H5R_DATASET_REGION = 1;

        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern int H5Rcreate(out RegRef newref, int loc_id, string name, int ref_type, int space_id);

        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern int H5Rdereference(int fileid, int ref_type, [In] ref RegRef refptr);

        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern int H5Rget_region(int fileid, int ref_type, [In] ref RegRef refptr);

        public static RegRef create(int loc_id, string name, int space_id)
        {
            RegRef r;

            CheckError(H5Rcreate(out r, loc_id, name, H5R_DATASET_REGION, space_id));

            return r;
        }

        /// <summary>
        /// Dereference a reference
        /// </summary>
        /// <param name="fileid"></param>
        /// <param name="refptr"></param>
        /// <returns>The dataset id pointed to by the reference</returns>
        public static int dereference(int fileid, RegRef refptr)
        {
            lock (HDFLock)
            {
                return CheckError(H5Rdereference(fileid, H5R_DATASET_REGION, ref refptr));
            }
        }

        /// <summary>
        /// Get the dataspace associated with a reference
        /// </summary>
        /// <param name="fileid"></param>
        /// <param name="refptr"></param>
        /// <returns></returns>
        public static int get_region(int fileid, RegRef refptr)
        {
            lock (HDFLock)
            {
                return CheckError(H5Rget_region(fileid, H5R_DATASET_REGION, ref refptr));
            }
        }
    }


    public class H5T : HDFGlue
    {
        /// <summary>
        /// C# doesn't support importing data from C DLLs, so we have to use the painful text_to_dtype operations
        /// </summary>
        public static readonly int NATIVE_INT = H5L.text_to_dtype("H5T_NATIVE_INT");

        public static readonly int NATIVE_SHORT = H5L.text_to_dtype("H5T_NATIVE_SHORT");
        public static readonly int NATIVE_USHORT = H5L.text_to_dtype("H5T_NATIVE_USHORT");
        public static readonly int NATIVE_DOUBLE = H5L.text_to_dtype("H5T_NATIVE_DOUBLE");

        /// <summary>
        /// The 'intel' byte ordering
        /// </summary>
        public static readonly int STD_I8LE = H5L.text_to_dtype("H5T_STD_I8LE");

        public static readonly int STD_I16LE = H5L.text_to_dtype("H5T_STD_I16LE");
        public static readonly int STD_I32LE = H5L.text_to_dtype("H5T_STD_I32LE");

        /// <summary>
        /// The "COMPOUND" type
        /// </summary>
        /*
        <datatype> ::= <atomic_type> | <compound_type> | <array_type> |
		<variable_length_type>

        <atomic_type> ::= <integer>  | <float>  | <time>      | <string> |
                          <bitfield> | <opaque> | <reference> | <enum>

        <integer> ::=  H5T_STD_I8BE     | H5T_STD_I8LE      |
                       H5T_STD_I16BE    | H5T_STD_I16LE     |
                       H5T_STD_I32BE    | H5T_STD_I32LE     |
                       H5T_STD_I64BE    | H5T_STD_I64LE     |
                       H5T_STD_U8BE     | H5T_STD_U8LE      |
                       H5T_STD_U16BE    | H5T_STD_U16LE     |
                       H5T_STD_U32BE    | H5T_STD_U32LE     |
                       H5T_STD_U64BE    | H5T_STD_U64LE     |
                       H5T_NATIVE_CHAR  | H5T_NATIVE_UCHAR  |
                       H5T_NATIVE_SHORT | H5T_NATIVE_USHORT |
                       H5T_NATIVE_INT   | H5T_NATIVE_UINT   |
                       H5T_NATIVE_LONG  | H5T_NATIVE_ULONG  |
                       H5T_NATIVE_LLONG | H5T_NATIVE_ULLONG

        <float> ::= H5T_IEEE_F32BE   | H5T_IEEE_F32LE     |
                    H5T_IEEE_F64BE   | H5T_IEEE_F64LE     |
                    H5T_NATIVE_FLOAT | H5T_NATIVE_DOUBLE  |
                    H5T_NATIVE_LDOUBLE

        <time> ::= TBD

        <string> ::= H5T_STRING { STRSIZE <strsize> ;
          STRPAD <strpad> ;
          CSET <cset> ;
          CTYPE <ctype> ;}
          
        <strsize> ::= <int_value> | H5T_VARIABLE
        <strpad> ::= H5T_STR_NULLTERM | H5T_STR_NULLPAD | H5T_STR_SPACEPAD
        <cset> ::= H5T_CSET_ASCII | H5T_CSET_UTF8
        <ctype> ::= H5T_C_S1 | H5T_FORTRAN_S1

        <bitfield> ::= TBD

        <opaque> ::= H5T_OPAQUE { OPQ_SIZE <opq_size>;
			           OPQ_TAG <opq_tag>; }
        opq_size ::= <int_value>
        opq_tag ::= “<string>”

        <reference> ::= Not supported

        <compound_type> ::= H5T_COMPOUND { <member_type_def>+ }
        <member_type_def> ::= <datatype> <field_name> <offset>opt ;
        <field_name> ::= “<identifier>”
        <offset> ::= : <int_value>

        <variable_length_type> ::= H5T_VLEN { <datatype> }

        <array_type> ::= H5T_ARRAY { <dim_sizes> <datatype> }
        <dim_sizes> ::= [<dimsize>] | [<dimsize>] <dim_sizes>
        <dimsize> ::= <int_value>

        <enum> ::= H5T_ENUM { <enum_base_type>; <enum_def>+ }
        <enum_base_type> ::= <integer>
        // Currently enums can only hold integer type data, but they may be 
        //expanded in the future to hold any datatype
        <enum_def> ::= <enum_symbol> <enum_val>;
        <enum_symbol> ::= “<identifier>”
        <enum_val> ::= <int_value>

        */
        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern int H5Tclose(int id);

        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern int H5Tget_size(int id);

        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern int H5Tset_size(int id, IntPtr size);

        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern int H5Tcopy(int id);


        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern int H5Tequal(int id1, int id2);

        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        public static extern int H5Tcreate(int h5type, int size);

        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        public static extern int H5Tinsert(int parent_id, string name, int offset, int member_id);

        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern int H5Tget_native_type(int dtype_id, int direction);

        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern int H5Tget_class(int dtype_id);

        public static void close(int id)
        {
            lock (HDFLock)
            {
                CheckError(H5Tclose(id));
            }
        }

        public static int get_size(int id)
        {
            lock (HDFLock)
            {
                return CheckError(H5Tget_size(id));
            }
        }

        public static int set_size(int id, int size)
        {
            lock (HDFLock)
            {
                IntPtr sz = new IntPtr(size);
                return CheckError(H5Tset_size(id, sz));
            }
        }


        public static bool equal(int id1, int id2)
        {
            lock (HDFLock)
            {
                return CheckError(H5Tequal(id1, id2)) > 0;
            }
        }

        public static int copy(int id)
        {
            lock (HDFLock)
            {
                return CheckError(H5Tcopy(id));
            }
        }

        public static int get_native_type(int id)
        {
            lock (HDFLock)
            {
                return CheckError(H5Tget_native_type(id, 0));
            }
        }

        public static int get_class(int id)
        {
            lock (HDFLock)
            {
                return CheckError(H5Tget_class(id));
            }
        }
    }

    public class H5Z : HDFGlue
    {
        public static int H5Z_FILTER_DEFLATE = 1;


        [DllImport(hdfPInvokeName, CallingConvention = CallingConvention.Cdecl)]
        private static extern int H5Zfilter_avail(int id);

        public static int filter_avail(int id)
        {
            lock (HDFLock)
                return CheckError(H5Zfilter_avail(id));
        }
    }
}

// ReSharper restore InconsistentNaming