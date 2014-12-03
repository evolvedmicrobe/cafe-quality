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
using System.Linq;
using System.Runtime.InteropServices;
using PacBio.Utils;

namespace PacBio.HDF
{

    /// <summary>
    /// Class containing various 'unsupported' interactions with HDF5. When you need to do something
    /// specia, possibly involving strings
    /// </summary>
    public static class HDFUtils
    {
        /// <summary>
        /// Special support for writing arrays of strings to HDF5 attributes.  Used for writing 'Content' attributes
        /// for pulse / base files
        /// </summary>
        /// <param name="target">Group / Dataset to add attribute to</param>
        /// <param name="name">Name of new attribute</param>
        /// <param name="strings">Array of strings to add</param>
        /// <returns>The newly created IAttribute</returns>
        public static IDataContainer WriteStringArrayAttribute(this IAttributeTarget target, string name,  Array strings)
        {
            if (target.GetAttribute(name) != null)
                target.DeleteAttribute(name);

            var dims = strings.Dimensions();

            var maxLen = 0;

            foreach(object o in strings)
            {
                string s = o as string;

                if (s != null && s.Length > maxLen)
                    maxLen = s.Length;
            }
            maxLen++;

            var t = HDFGlue.GetHDFGlobal("H5T_C_S1_g");
            t = H5T.copy(t);
            var typ = new HDFDatatype(t);

            H5T.set_size(typ.Id, maxLen);

            var dataspace = new HDFDataspace(target.File, dims, dims);
            var attr = target.CreateAttribute(name, typ, dataspace) as HDFAttribute;
            
            if(attr == null)
                throw new ApplicationException("Unable to create HDF5 attribute");

            // Allocate memory for strings
            var buf = Marshal.AllocHGlobal((int)(maxLen*dims.Aggregate((a, b) => a*b)));
            var bp = buf;

            
            foreach (var l in strings.EnumerateIndicies())
            {
                var source = ((string)strings.GetValue(l)).ToCharArray().Map(c => (byte)c);
                var byteBuf = new byte[maxLen];
                Array.Copy(source, byteBuf, source.Length);

                Marshal.Copy(byteBuf, 0, bp, byteBuf.Length);
                bp = bp.Offset(maxLen);
            }

            H5A.write(attr.Id, typ.Id, buf);
            Marshal.FreeHGlobal(buf);
            return attr;
        }


        /// <summary>
        /// Add a byte count offset to an IntPtr
        /// </summary>
        /// <param name="src"></param>
        /// <param name="offset"></param>
        /// <returns></returns>
        public static IntPtr Offset(this IntPtr src, int offset)
        {
            switch (IntPtr.Size)
            {
                case 4:
                    return new IntPtr(src.ToInt32() + offset);
                case 8:
                    return new IntPtr(src.ToInt64() + offset);
                default:
                    throw new NotSupportedException("Surprise!  This is running on a machine where pointers are " + IntPtr.Size + " bytes and arithmetic doesn't work in C# on them.");
            }
        }

        /// <summary>
        /// From a group node, find all elements of type desired,
        /// return a list of elements. 
        /// </summary>
        /// <param name="root"></param>
        /// <param name="type"></param>
        /// <returns></returns>
        public static List<string> FindAllElements(IGroup root, Type type)
        {
            var nodeList = new List<string>();
            
            // Construct the node list for elements of input type
            if (root != null)
            {
                var nodes = root.GetChildren();
                nodes.Where(d => (d.GetType() == type)).ForEach(v => nodeList.Add(v.Name));
            }

            return nodeList;
        }
    }
}
