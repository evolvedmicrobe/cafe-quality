using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using PacBio.HDF;
using PacBio.Utils;

namespace PacBio.IO.ICD
{
    /// <summary>
    /// An entry in an Interface Control Document for and HDF file.
    /// </summary>
    public class HdfIcdEntry : Attribute
    {
        /// <summary>
        /// A description of the object
        /// </summary>
        public string Detail { get; set; }

        /// <summary>
        /// The full HDF path to the object
        /// </summary>
        public string Path { get; set; }

        /// <summary>
        /// Units or Encoding (where relevant) of the data
        /// </summary>
        public object Units { get; set; }

        /// <summary>
        /// Flag to indicate that this entry was found in the ICD,
        /// as opposed to inferred from the file
        /// </summary>
        public bool InIcdDoc { get; set; }

        /// <summary>
        /// Flag to indicate a deprecated field, not supported by PacBio RS.
        /// </summary>
        public bool Deprecated { get; set; }
    }

    /// <summary>
    /// An Interface Control Document as a Dictionary of named HdfIcdEntry elements.
    /// </summary>
    public class Icd : Dictionary<string, HdfIcdEntry>
    {
        /// <summary>
        /// Concatenate a partial ICD with this one.
        /// </summary>
        /// <param name="part"></param>
        public void Concat(Icd part)
        {
            foreach (var entry in part)
            {
                if (!ContainsKey(entry.Key))
                    Add(entry.Key, entry.Value);
            }
        }
    }

    /// <summary>
    /// A Hierarchical Data Format ICD, based on a Reader type T.
    /// The reader class should be annotated with an HdfIcdEntry for the containing group, and
    /// it can contain HDF dataset or HDF attribute entries as HdfIcdEntry attributes on public
    /// members that directly access the corresponding file elements.  Alternatively, "manual"
    /// ICD entries can be made in the Icd class contained in Reader type T.  These should be
    /// returned in a GetIcd() method, to be invoked through reflection. The design pattern is:
    /// 
    /// class Reader
    /// {
    ///     class Icd : HdfIcd{Reader}
    ///     {
    ///         // ...
    ///     }
    /// 
    ///     public static Data.Icd GetIcd() { return new Icd( /* flatten = */ true); }
    /// }
    /// </summary>
    /// <typeparam name="T"></typeparam>
    public class HdfIcd<T> : Icd
    {
        // Flatten the path hierarchy to name-only in the dictionary
        private readonly bool flatten;

        // Entry of the containing group, which should be made on the Reader (T) class
        private readonly HdfIcdEntry groupEntry;

        protected HdfIcd(IEnumerable<HdfIcdEntry> manEntries = null, bool flatten = false)
        {
            this.flatten = flatten;

            // Check the class for an annotation
            var readerType = typeof(T);

            var cattr = readerType.GetCustomAttributes(typeof(HdfIcdEntry), true);
            if (cattr.Length == 1)
            {
                groupEntry = (HdfIcdEntry)cattr[0];
            }
            else
            {
                // Require that the designated type T is annotated with the corresponding group info
                throw new ApplicationException("Class must contain the HdfIcdEntry attribute for the group");
            }

            // Manual entries
            if (manEntries != null)
                manEntries.ForEach(Append);

            // Member entries
            AppendMemberEntries(readerType);

            // Add the reader group as an entry
            Append(groupEntry);
        }

        private void AppendMemberEntries(Type type)
        {
            var members = type.GetMembers();
            foreach (var mi in members)
            {
                var attr = mi.GetCustomAttributes(typeof(HdfIcdEntry), true);
                if (attr.Length == 0) continue;

                var icdEntry = (HdfIcdEntry)attr[0];

                // If the Location is not specified, use the member name
                if (icdEntry.Path == null)
                {
                    icdEntry.Path = mi.Name;
                }
                else if (!icdEntry.Path.StartsWith("/") && icdEntry.Path.EndsWith("/"))
                {
                    // Alternatively, a relative path can be provided.
                    icdEntry.Path = Path.Combine(icdEntry.Path, mi.Name);
                }

                Append(icdEntry);
            }
        }

        private void Append(HdfIcdEntry entry)
        {
            // Re-define the Location to be the full path
            if (!entry.Path.StartsWith("/"))
                entry.Path = Path.Combine(groupEntry.Path, entry.Path);

            var name = flatten ? Path.GetFileName(entry.Path) : entry.Path;

            if (!ContainsKey(name))
                Add(name, entry);
        }
    }
}