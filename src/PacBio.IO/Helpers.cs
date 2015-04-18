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
using System.Globalization;
using System.IO;
using System.Linq;
using PacBio.HDF;

namespace PacBio.IO
{
    /// <summary>
    /// Some helpful extension methods
    /// </summary>
    public static class Helpers
    {
        /// <summary>
        /// Open the URI and check if it is an HDF multi-part index file.
        /// If it is, return the URIs of the part files; otherwise return null.
        /// </summary>
        /// <param name="filename"></param>
        /// <returns></returns>
        public static IEnumerable<string> GetMultiPartUris(string filename)
        {
            if (!System.IO.File.Exists(filename))
            {
                throw new FileNotFoundException(String.Format("Multipart HDF not found: {0}",
                                                              filename), filename);
            }

            var path = Path.GetDirectoryName(filename);

            using (var chunks = HDFFile.Open(filename, FileMode.Open, FileAccess.Read))
            {
                var mpg = (HDFGroup) chunks.GetChild("MultiPart");
                if (mpg != null)
                {
                    // Now it must be a valid multi-part file
                    var parts = (IDataset) mpg.GetChild("Parts");
                    var names = (string[]) parts.Read();

                    return names.Select(n => Path.Combine(path, n));
                }

                return null;
            }
        }

        public static DateTime GetDateAttribute(IAttributeTarget hdfGroup, string attributeName)
        {
            var dateAttribute = hdfGroup.GetAttribute("DateCreated");

            // If we can't find the attribute, get the file date
            if (dateAttribute == null)
            {
                try
                {
                    return FileTime(hdfGroup);

                }
                catch (Exception)
                {
                    return new DateTime(1900, 1, 1);
                }
            }

            DateTime dateCreated;
            var dateString = dateAttribute.ReadSingleton<string>();
            var parsed = DateTime.TryParse(dateString, null, DateTimeStyles.RoundtripKind, out dateCreated);

            // Try it this way for legacy Matlab produced data.
            if(!parsed)
            {
                parsed = DateTime.TryParseExact(dateString, "yyyyMMdd'T'HHmmss", null, DateTimeStyles.None, out dateCreated);
            }

            if (parsed)
                return dateCreated;

            return new DateTime(1900, 1, 1);
        }


        static DateTime FileTime(IChunkElement hdfGroup)
        {
            var fileName = hdfGroup.File.Name;
            return File.GetCreationTime(fileName);
        }
    }
}
