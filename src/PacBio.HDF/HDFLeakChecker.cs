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
using System.Text;

namespace PacBio.HDF
{
    /// <summary>
    /// Check for leaks in HDF (by creating a snapshot with a count of open object IDs)
    /// </summary>
    /// When we get disposed we will check for leaks (useful with the using construct)
    public class HDFLeakChecker : IDisposable
    {
        HDFFile src;

        private H5F.ObjTypes[] tocheck = {
                                             H5F.ObjTypes.OBJ_ALL,
                                             H5F.ObjTypes.OBJ_ATTR, 
                                             H5F.ObjTypes.OBJ_DATASET, 
                                             H5F.ObjTypes.OBJ_DATATYPE,
                                             H5F.ObjTypes.OBJ_GROUP
                                         }; 
        // Don't check files for now.
        // H5F.ObjTypes.OBJ_FILE };

        IDictionary<H5F.ObjTypes, int> data = new Dictionary<H5F.ObjTypes, int>();

        /// <summary>
        /// We are not multithread safe yet
        /// </summary>
        static bool Enabled = false;

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="src">the file to watch</param>
        public HDFLeakChecker(HDFFile src) : this(src, false)
        {
        }

        /// <summary>
        /// Setup and HDF5 leak checker
        /// </summary>
        /// <param name="src">The file to watch</param>
        /// <param name="forceChecking">Override the static Enabled flag and get object counts in constructor</param>
        public HDFLeakChecker(HDFFile src, bool forceChecking)
        {
            if (Enabled || forceChecking)
            {
                this.src = src;

                if (src.Id < 0)
                    return;

                foreach (var typ in tocheck)
                    data.Add(typ, H5F.get_obj_count(src.Id, typ));
            }
        }

        /// <summary>
        /// Compare the current object count to the original snapshot
        /// </summary>
        public void CheckLeaks()
        {
            var newer = new HDFLeakChecker(src);

            foreach (var pair in data)
            {
                var newcount = newer.data[pair.Key];
                if (pair.Value < newcount)
                {
                    throw new ApplicationException(string.Format("HDF leak of {0}, old={1}, new={2}", pair.Key, pair.Value, newcount));
                }
            }

            // Coverity complains about this because this class implements IDisposable
            // HOWEVER, it uses the disposable pattern to finish leak detection and then
            // calls THIS routine to perform the dispose. So disposing here is probably
            // a bad idea.
        }

        public override string ToString()
        {
            StringBuilder b = new StringBuilder();

            b.Append("HDF Object counts: ");

            foreach (var tuple in data)
                b.AppendFormat("{0}={1} ", tuple.Key, tuple.Value);

            return b.ToString();
        }

        #region IDisposable Members

        public void Dispose()
        {
            Dispose(true);
            GC.SuppressFinalize(this);
        }

        protected virtual void Dispose(bool disposing)
        {
            if (Enabled)
                CheckLeaks();
        }

        #endregion
    }
}
