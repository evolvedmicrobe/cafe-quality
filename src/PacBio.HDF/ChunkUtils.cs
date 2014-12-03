using System;
using PacBio.Utils;

namespace PacBio.HDF
{
    /// <summary>
    /// High level API that works with any chunk provider
    /// </summary>
    public static class ChunkUtils
    {
        /// <summary>
        /// Given a root node and a path, creat any needed groups and return the modified path
        /// </summary>
        /// <param name="root"></param>
        /// <param name="path"></param>
        /// <returns></returns>
        static public IGroup ExpandPath(this IGroup root, ref string path)
        {
            // If the user starts their search at /, then we want to handle that case specially
            bool startAtRoot = false;
            if(path[0] == '/')
            {
                path = path.Substring(1);
                startAtRoot = true;
            }

            string[] dirs = path.Split(new char[] {'/'});

            // Start our search at the root of the hierarchy if needed

            if (startAtRoot && dirs.Length > 0)
                dirs[0] = "/" + dirs[0];

            int numsubdirs = dirs.Length - 1;   // We create a directory for every node about the _last_ node

            int cursubdir = 0;
            while(numsubdirs-- > 0)
            {
                string dirname = dirs[cursubdir++];

                IGroup child = root.CreateGroup(dirname);

                root = child;   // FIXME - we are leaking memory here with HDF (need to dispose these groups as they are used)
            }

            // The last part of the path is the name of our new element
            path = dirs[cursubdir];

            return root;
        }


        /// <summary>
        /// Create or open a group as needed
        /// </summary>
        /// <param name="parent"></param>
        /// <param name="path"></param>
        /// <returns></returns>
        public static IGroup Create(this IGroup parent, string path)
        {
            using(IGroup mygroup = ExpandPath(parent, ref path))  // Coverity RESOURCE_LEAK
                return mygroup.CreateGroup(path);
        }


        /// <summary>
        /// Create a dataspace with the correct dimensions to hold the specified type
        /// </summary>
        /// <param name="parent"></param>
        /// <param name="obj"></param>
        /// <returns></returns>
        public static IDataspace CreateDataspace(this IGroup parent, object obj)
        {
            Array arr;

            if ((arr = obj as Array) != null)
            {
                long [] dims = new long[arr.Rank];

                for (int i = 0; i < dims.Length; i++)
                    dims[i] = arr.GetLength(i);

                // We might need to swap the array dims
                // if (!parent.File.RowMajorOrder)
                //    CollectionTools.Swap(ref dims);

                return parent.File.CreateDataspace(dims, dims);
            }
            else
                // Create a scalar space
                return parent.File.CreateDataspace();
        }

        /// <summary>
        /// Simply write an object to a chunk file - creating datasets etc... as needed
        /// </summary>
        /// <param name="parent"></param>
        /// <param name="path"></param>
        /// <param name="data"></param>
        public static void Create(this IGroup parent, string path, object data)
        {
            IGroup mygroup = ExpandPath(parent, ref path);

            // First see if we can reuse an existing dataset
            IDataset dset = (IDataset)mygroup.GetChild(path);

            // Guess we need to make it?
            if (dset == null)
                using (IDatatype dtype = parent.File.CreateDatatype(data.GetType()))
                using (var dspace = CreateDataspace(parent, data))
                    dset = mygroup.CreateDataset(path, dtype, dspace);

            using (dset)
                dset.Write(data);

        }
        

        /// <summary>
        /// Add a simply typed attribute to the IGroup
        /// </summary>
        /// <param name="target">Object to add attribute to</param>
        /// <param name="name">Name of the attribute</param>
        /// <param name="o">Attribute data. Must be a simple datatype or a string.</param>
        public static IDataContainer InsertAttribute(this IAttributeTarget target, string name, object o)
        {
            if (target.GetAttribute(name) != null)
                target.DeleteAttribute(name);


            var t = o.GetType();
            var hdfType = target.File.CreateDatatype(t);

            IDataspace space;

            if (o is Array)
            {
                var a = o as Array;
                space = target.File.CreateDataspace(a.Dimensions(), a.Rank.Fill(-1L));
            }
            else
            {
                space = target.File.CreateDataspace();
            }

            var attr = target.CreateAttribute(name, hdfType, space);
            attr.Write(o);
            return attr;
        }




        /// <summary>
        /// Add a simply typed attribute to the IGroup
        /// </summary>
        /// <param name="group">IGroup to add the attribute to</param>
        /// <param name="name">Name of the attribute</param>
        /// <param name="o">Attribute data. Must be a simple an array of a simple datatype</param>
        public static IDataset AddDataset(this IGroup group, string name, object o)
        {
            if (group.GetChild(name) != null)
                group.RemoveChild(name);

            var t = o.GetType();
            var a = o as Array;
            
            if (a == null)
                throw new ArgumentException("Must pass an array of a simple datatype");

            var space = group.File.CreateDataspace(a.Dimensions(), a.Rank.Fill(-1L));

            var ds = group.CreateDataset(name, group.File.CreateDatatype(t), space);
            ds.Write(o);
            return ds;
        }


        /// <summary>
        /// Deal with singletons that may have either been stored as a single value or in a 1 element array
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="data"></param>
        /// <returns></returns>
        public static T ReadSingleton<T>(this IDataContainer data)
        {
            var o = data.Read();

            if (o is T)
                return (T) o;

            if (o is T[])
                return ((T[]) o)[0];

            throw new ArgumentException(
                String.Format("Stored type is {0}, requested type is {1}", o.GetType(), typeof(T)));
        }

        public static object ReadDataset(this IGroup group, string datasetName)
        {
            var ds = (IDataset)group.GetChild(datasetName);

            if (ds == null)
                return null;        // oops - nothing here by that name

            return ds.Read();
        }

    }

}

