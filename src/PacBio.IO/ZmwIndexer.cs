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
using PacBio.HDF;
using PacBio.Utils;

namespace PacBio.IO
{
    public struct XYPoint
    {
        public int X;
        public int Y;
    }

    /// <summary>
    /// Manages the mappings between Zmw HoleNumber, Zmw X/Y, and the raw index of a zmw in an HDF5 file
    /// </summary>
    public class ZmwIndexer
    {
        /// <summary>
        /// The HoleNumbers of the Zmws represented in this file
        /// </summary>
        public int[] HoleNumber { get; private set; }

        /// <summary>
        /// The look vector of the data
        /// </summary>
        public Int16[] HoleChipLook { get; private set; }

        /// <summary>
        /// The X/Y coordinates of the Zmws represented in this file
        /// </summary>
        public XYPoint[] HoleXY { get; private set; }

        /// <summary>
        /// The ZmwType of the Zmws stored in this file
        /// </summary>
        public ZmwType[] ZmwType { get; private set; }

        internal Dictionary<int, int> HoleNumberMap;
        internal Dictionary<int, Dictionary<int, int>> xyMap
        {
            get
            {
                if (xyMap == null)
                    InitXYMap();

                return _xyMap;
            }
        }

        private Dictionary<int, Dictionary<int, int>> _xyMap; 
        
        /// <summary>
        /// Number of Zmw stored in this file
        /// </summary>
        internal int NumZmws { get; private set; }
        
        public ZmwIndexer(int[] holeNumber, XYPoint[] holeXY)
        {
            HoleXY = holeXY;
            NumZmws = HoleXY.GetLength(0);

            HoleNumber = holeNumber;
            ZmwType = holeNumber.Length.Fill(i => IO.ZmwType.Sequencing);

            InitIndexes();
        }

        /// <summary>
        /// Set up a ZmwIndexer based on a ZMW group
        /// </summary>
        /// <param name="zmwGroup"></param>
        public ZmwIndexer(IGroup zmwGroup)
        {
            var holexy = (Int16[,])((IDataset)zmwGroup.GetChild("HoleXY")).Read();
            NumZmws = holexy.GetLength(0);
            HoleXY = NumZmws.Fill(h => new XYPoint { X = holexy[h, 0], Y = holexy[h, 1] }).ToArray();

            var hn = (uint[])((IDataset)zmwGroup.GetChild("HoleNumber")).Read();
            HoleNumber = hn.Map(i => (int)i);

            // Check for a ZMW look field -- this will only exist on 150k files
            // If it doesn't exist, just fill it out with zeros.
            var lookDataset = (IDataset) zmwGroup.GetChild("HoleLook");
            if (lookDataset != null)
            {
                HoleChipLook = (Int16[]) lookDataset.Read();
            }
            else
            {
                HoleChipLook = NumZmws.Fill(i => (Int16)0);
            }

            var zt = ((IDataset)zmwGroup.GetChild("HoleStatus"));
            if (zt != null)
            {
                // First update the ZmwType array with the contents of the status field
                var zmwTypeLookup = zt.GetAttribute("LookupTable");

                if (zmwTypeLookup != null)
                    IO.ZmwType.SetupTypes((string[])zmwTypeLookup.Read());

                ZmwType = ((byte[])zt.Read()).Map(i => (ZmwType)i);
            }
            else
            {
                ZmwType = NumZmws.Fill(i => IO.ZmwType.Sequencing);
            }


            InitIndexes();
        }

        private void InitIndexes()
        {
            // Map from HoleNumber to Zmw index
            HoleNumberMap = HoleNumber.Select((hn, i) => new { Key = hn, Value = i }).ToDictionary(e => e.Key,
                                                                                                 e => e.Value);
        }

        private void InitXYMap()
        {            
            // Map from X/Y to Zmw index
            _xyMap = new Dictionary<int, Dictionary<int, int>>();

            HoleXY.ForEach((idx, xy) =>
            {
                Dictionary<int, int> innerMap;
                var s = xyMap.TryGetValue(xy.X, out innerMap);

                if (!s)
                    xyMap[xy.X] = new Dictionary<int, int>();

                xyMap[xy.X][xy.Y] = idx;
            });
        }

        public int GetIndexByHoleNumber(int holeNumber)
        {
            return HoleNumberMap[holeNumber];
        }

        public int GetIndexByHoleXY(int x, int y)
        {
            return xyMap[x][y];
        }
    }
}