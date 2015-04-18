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
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using PacBio.HDF;
using PacBio.Utils;

namespace PacBio.IO
{
    /// <summary>
    /// A type to describe a generic 1D array slice, consistent with the HDF5
    /// hyperslab conventions: Stride is the number of elements between Block
    /// start indices; Count corresponds to the number of blocks (the number of
    /// elements is Block * Count).
    /// </summary>
    public class HdfSlice
    {
        /// <summary>
        /// The array index of the first block
        /// </summary>
        public int Start { get; set; }

        /// <summary>
        /// The number of elements between block start indices
        /// </summary>
        public int Stride { get; set; }

        /// <summary>
        /// The block size
        /// </summary>
        public int Block { get; set; }

        /// <summary>
        /// The number of blocks
        /// </summary>
        public int Count { get; set; }

        public HdfSlice()
        {
        }

        /// <summary>
        /// Copy constructor
        /// </summary>
        /// <param name="x"></param>
        public HdfSlice(HdfSlice x)
        {
            Start = x.Start;
            Stride = x.Stride;
            Block = x.Block;
            Count = x.Count;
        }
    }

    /// <summary>
    /// An interface for ZMW range selection by strided blocks of HoleNumber.
    /// </summary>
    public interface IZmwRange
    {
        /// <summary>
        /// The first (inclusive) HoleNumber.
        /// </summary>
        int Start { get; set; }

        /// <summary>
        /// The number of blocks to include.
        /// </summary>
        int Count { get; set; }

        /// <summary>
        /// The block size.
        /// </summary>
        int Block { get; set; }

        /// <summary>
        /// HoleNumber distance between the start of block i and the start of block i+1.
        /// </summary>
        int Stride { get; set; }
    }

    /// <summary>
    /// A type to denote a 1D slice over a HoleNumber range. When applied to a list
    /// corresponding to non-consecutive HoleNumbers, the returned items should include
    /// only those elements that overlap with the ordinary HdfSlice on the list of
    /// consecutive HoleNumbers.
    /// </summary>
    public class ZmwRange : HdfSlice, IZmwRange
    {
        /// <summary>
        /// Default constructor to select all available HoleNumbers in a file.
        /// </summary>
        public ZmwRange()
        {
            Start = 0;
            Stride = 1;
            Count = -1;
            Block = 1;
        }

        public override string ToString()
        {
            return String.Format("Start={0}, Stride={1}, Count={2}, Block={3}, ", Start, Stride, Count, Block);
        }
    }

    /// <summary>
    /// An abstract class to provide implementors of IDataSource the basic logic for accessing
    /// items, esp. using a ZmwRange. Implementations must set the ZmwSource and provide,
    /// minimally, the List indexing method [].  Override other generic interface methods as
    /// required.
    /// </summary>
    /// <typeparam name="T"></typeparam>
    public abstract class DataSource<T> : IDataSource<T>
    {
        // Implementation must set the ZmwSource
        public virtual IZmwSource ZmwSource { get; protected set; }

        // The first (lowest) HoleNumber provided by the source
        public virtual int FirstHoleNumber { get { return ZmwSource[0].HoleNumber; } }

        // The last (highest) HoleNumber provided by the source
        public virtual int LastHoleNumber { get { return ZmwSource[Count - 1].HoleNumber; } }

        #region ISourceIdentifier

        public virtual string FileName { get; protected set; }

        public virtual DateTime DateCreated
        {
            get { throw new NotImplementedException(); }
        }

        public virtual string SoftwareVersion
        {
            get { return "Unknown Version"; }
        }

        public virtual string ChangelistID
        {
            get { return "Unknown Revision"; }
        }

        #endregion

        #region IDataSource<T>

        public IMovieMetadata Movie
        {
            get { return ZmwSource; }
        }

        public virtual T ByHoleNumber(int holeNum)
        {
            return this[ZmwSource.GetIndexByHoleNumber(holeNum)];
        }

        public virtual T ByXY(int x, int y)
        {
            return this[ZmwSource.GetIndexByHoleXY(x, y)];
        }

        /// <summary>
        /// Access data objects by hole-number range.
        /// </summary>
        /// <param name="range">Range parameter. Count=-1 means all.</param>
        /// <returns>Enumerable list of data. Range is Start inclusive.</returns>
        public IEnumerable<T> ByHoleNumberRange(IZmwRange range)
        {
            // Null range will select all.
            range = range ?? new ZmwRange();

            // Requires increasing holeNumber with index
            int lastHoleNumber = ZmwSource[Count - 1].HoleNumber;

            int lower = range.Start;
            int upper = (range.Count <= 0 ? lastHoleNumber + 1 : range.Start + range.Stride * range.Count);

            if (upper == lower)
                throw new ArgumentException("Invalid range: " + range);

            // Keep it readable
            // ReSharper disable LoopCanBeConvertedToQuery
            foreach (var zmw in ZmwSource)
            {
                int hn = zmw.HoleNumber;
                int blockEnd = Math.Min(upper, range.Start + ((hn - range.Start) / range.Stride) * range.Stride + range.Block);

                if (lower <= hn && hn < blockEnd)
                    yield return ByHoleNumber(hn);
            }
            // ReSharper restore LoopCanBeConvertedToQuery
        }

        public virtual void ValidateSource() { throw new NotImplementedException(); }

        #endregion

        #region IEnumerable

        public virtual IEnumerator<T> GetEnumerator()
        {
            return ZmwSource.Select(z => ByHoleNumber(z.HoleNumber)).GetEnumerator();
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            return GetEnumerator();
        }

        #endregion

        #region ICollection

        public virtual void Add(T item)
        {
            throw new NotImplementedException();
        }

        public virtual void Clear()
        {
            throw new NotImplementedException();
        }

        public virtual bool Contains(T item)
        {
            throw new NotImplementedException();
        }

        public virtual void CopyTo(T[] array, int arrayIndex)
        {
            throw new NotImplementedException();
        }

        public virtual bool Remove(T item)
        {
            throw new NotImplementedException();
        }

        public virtual int Count
        {
            get { return ZmwSource.ZmwIndexer.NumZmws; }
        }

        public bool IsReadOnly
        {
            get { return true; }
        }

        #endregion

        #region IList

        public virtual int IndexOf(T item)
        {
            throw new NotImplementedException();
        }

        public virtual void Insert(int index, T item)
        {
            throw new NotImplementedException();
        }

        public virtual void RemoveAt(int index)
        {
            throw new NotImplementedException();
        }

        public abstract T this[int index] { get; set; }

        #endregion

        #region IDisposable

        protected abstract void Dispose(bool disposing);

        public void Dispose()
        {
            Dispose(true);
            GC.SuppressFinalize(this);
        }

        ~DataSource()
        {
            Dispose(false);
        }

        #endregion
    }



    public class ZmwMultiPartSource : MultiPartSource<ISequencingZmw>, IZmwSource
    {
        private ZmwIndexer indexer { get; set; }

        public ZmwMultiPartSource(IEnumerable<IZmwSource> zmwSources) :
            base(zmwSources.Select(v => v as DataSource<ISequencingZmw>))
        {
            var holeNumber = parts.SelectMany(p => p.ZmwSource.ZmwIndexer.HoleNumber).ToArray();
            var holeXY = parts.SelectMany(p => p.ZmwSource.ZmwIndexer.HoleXY).ToArray();
            indexer = new ZmwIndexer(holeNumber, holeXY);
        }

        public string MovieName
        {
            get { return parts[0].Movie.MovieName; }
        }

        public string RunCode
        {
            get { return parts[0].Movie.RunCode; }
        }

        public string InstrumentName
        {
            get { return parts[0].Movie.InstrumentName; }
        }

        public uint InstrumentId
        {
            get { return parts[0].Movie.InstrumentId; }
        }

        public uint PlatformId
        {
            get { return parts[0].Movie.PlatformId; }
        }

        public string PlatformName
        {
            get { return parts[0].Movie.PlatformName; }
        }
        
        public float FrameRate
        {
            get { return parts[0].Movie.FrameRate; }
        }

        public uint NumFrames
        {
            get { return parts[0].Movie.NumFrames; }
        }

        public int LaserOnFrame
        {
            get { return parts[0].Movie.LaserOnFrame; }
        }

        public int HotStartFrame
        {
            get { return parts[0].Movie.HotStartFrame; }
        }

        public char[] BaseMap
        {
            get { return parts[0].Movie.BaseMap; }
        }

        public IGroup ScanDataGroup
        {
            get { return parts[0].Movie.ScanDataGroup; }
        }

        public string BindingKit
        {
            get { return parts[0].Movie.BindingKit; }
        }

        public string SequencingKit
        {
            get { return parts[0].Movie.SequencingKit; }
        }

        public string BaseCallerChangelistID
        {
            get { return parts[0].Movie.BaseCallerChangelistID; }
        }

        public PacBio.Data.ChemistryBarcodeTriple ChemistryBarcode
        {
            get { return parts[0].Movie.ChemistryBarcode; }
        }

        public string SequencingChemistry
        {
            get { return parts[0].Movie.SequencingChemistry; }
        }

        public int GetIndexByHoleNumber(int holeNumber)
        {
            return HoleNumToIndex(holeNumber);
        }

        public int GetIndexByHoleXY(int x, int y)
        {
            return XYToIndex(x, y);
        }

        public ZmwIndexer ZmwIndexer
        {
            get { return indexer; }
        }
    }

    /// <summary>
    /// A base class for constructing multi-part data sources, i.e. readers 
    /// (e.g. pls.h5, bas.h5, etc), where the source consists of mulitple parts.
    /// </summary>
    /// <typeparam name="TOutputType"></typeparam>
    public class MultiPartSource<TOutputType> : IDataSource<TOutputType> where TOutputType : class
    {
        #region Members

        internal readonly DataSource<TOutputType>[] parts;
        private readonly ZmwIndexer[] indexers;

        private readonly int numZmws;

        /// <summary>
        /// The list of URIs as ordered, by ascending hole-number chunks.
        /// </summary>
        public readonly string[] OrderedFileNames;

        // Provide access to a corresponding IZmwSource
        public IZmwSource ZmwSource { get; private set; }

        #endregion

        #region Structors

        protected MultiPartSource(IEnumerable<DataSource<TOutputType>> sources)
        {
            // Store the source parts, ordered by starting HoleNumber
            parts = sources.OrderBy(p => p.FirstHoleNumber).ToArray();
            
            // Provide access to the URIs of the parts in order
            OrderedFileNames = parts.Select(s => s.FileName).ToArray();

            indexers = parts.Length.Fill(i => parts[i].ZmwSource.ZmwIndexer);

            // Set the total number of ZMWs
            foreach (var p in parts)
                numZmws += p.Count;

            // Verify that the ordering contract is met. This check assumes (but does not verify)
            // that within the parts, the hole numbers are increasing.
            for (int i = 1; i < parts.Length; i++)
            {
                if (parts[i].FirstHoleNumber <= parts[i - 1].LastHoleNumber)
                    throw new ApplicationException("Input source parts are not valid HoleNumber ranges.");
            }

            // Construct the ZmwSource
            if (parts[0] as IZmwSource == null)
            {
                ZmwSource = new ZmwMultiPartSource(parts.Select(p => p.ZmwSource));
            }
            else
            {
                ZmwSource = (IZmwSource) this;
            }
        }

        private bool disposed = false;

        protected virtual void Dispose(bool disposing)
        {
            if (disposing && !disposed)
            {
                parts.ForEach(p => p.Dispose());
                disposed = true;
            }
        }

        public void Dispose()
        {
            Dispose(true);
            GC.SuppressFinalize(this);
        }

        ~MultiPartSource()
        {
            Dispose(false);
        }

        #endregion

        /// <summary>
        /// Return the array of hole numbers contained in the given part.
        /// </summary>
        /// <param name="partIndex"></param>
        /// <returns></returns>
        public int[] GetHoleNumbersByPart(int partIndex)
        {
            return indexers[partIndex].HoleNumber;
        }

        #region ISourceIdentifier

        /// <summary>
        /// Returns the URI of the first part.
        /// </summary>
        public string FileName
        {
            get { return parts[0].FileName; }
        }

        /// <summary>
        /// Returns creation date of the first part.
        /// </summary>
        public DateTime DateCreated
        {
            get { return parts[0].DateCreated; }
        }

        public string SoftwareVersion
        {
            get { return parts[0].SoftwareVersion; }
        }

        public string ChangelistID
        {
            get { return parts[0].ChangelistID; }
        }

        #endregion

        #region IDataSource<T>

        public IMovieMetadata Movie
        {
            get { return parts[0].Movie; }
        }

        public TOutputType ByHoleNumber(int holeNum)
        {
            for (int i = 0; i < indexers.Length; i++)
            {
                var idx = indexers[i];

                if (idx.HoleNumberMap.ContainsKey(holeNum))
                    return parts[i].ByHoleNumber(holeNum);
            }

            throw new ApplicationException("Hole number not found");
        }

        public TOutputType ByXY(int x, int y)
        {
            for (int i = 0; i < indexers.Length; i++)
            {
                var idx = indexers[i];

                if (idx.xyMap.ContainsKey(x))
                {
                    var yMap = idx.xyMap[x];
                    if (yMap.ContainsKey(y))
                        return parts[i].ByXY(x, y);
                }
            }

            throw new ApplicationException("Hole [X,Y] coordinates not found");
        }

        public IEnumerable<TOutputType> ByHoleNumberRange(IZmwRange range)
        {
            // Could be null, that's ok.
            range = range ?? new ZmwRange();

            // Rely on the fact that parts are ordered on construction;
            // HoleNumbers are in ascending order in all files.
            //
            int last = parts.Length - 1;
            int lastHoleNumber = parts[last].LastHoleNumber;

            int lower = range.Start;
            int upper = (range.Count <= 0 ? lastHoleNumber + 1 : range.Start + range.Stride * range.Count);

            if (upper == lower)
                throw new ArgumentException("Invalid range: " + range);

            // Keep it readable
            // ReSharper disable LoopCanBeConvertedToQuery
            foreach (var part in parts)
            {
                foreach (var zmw in part.ZmwSource)
                {
                    int hn = zmw.HoleNumber;
                    int blockEnd = Math.Min(upper, range.Start + (hn/range.Stride)*range.Stride + range.Block);

                    if (lower <= hn && hn < blockEnd)
                        yield return ByHoleNumber(hn);
                }
            }
            // ReSharper restore LoopCanBeConvertedToQuery
        }

        protected int HoleNumToIndex(int holeNum)
        {
            int offset = 0;
            foreach (var idx in indexers)
            {
                if (idx.HoleNumberMap.ContainsKey(holeNum))
                    return offset + idx.HoleNumberMap[holeNum];

                offset += idx.NumZmws;
            }

            throw new ApplicationException("Hole number not found");
        } 

        protected int XYToIndex(int x, int y)
        {
            int offset = 0;
            foreach (var idx in indexers)
            {
                if (idx.xyMap.ContainsKey(x))
                {
                    var yMap = idx.xyMap[x];
                    if (yMap.ContainsKey(y))
                        return offset + yMap[y];
                }

                offset += idx.NumZmws;
            }

            throw new ApplicationException("Hole [X,Y] coordinates not found");
        }

        public virtual void ValidateSource() { throw new NotImplementedException(); }

        #endregion

        #region IEnumerable

        public IEnumerator<TOutputType> GetEnumerator()
        {
            return ByHoleNumberRange(null).GetEnumerator();
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            return GetEnumerator();
        }

        #endregion

        #region ICollection

        public void Add(TOutputType item)
        {
            throw new NotImplementedException();
        }

        public void Clear()
        {
            throw new NotImplementedException();
        }

        public bool Contains(TOutputType item)
        {
            throw new NotImplementedException();
        }

        public void CopyTo(TOutputType[] array, int arrayIndex)
        {
            var evArr = ByHoleNumberRange(null).ToArray();
            evArr.CopyTo(array, arrayIndex);
        }

        public bool Remove(TOutputType item)
        {
            throw new NotImplementedException();
        }

        public int Count
        {
            get { return numZmws; }
        }

        public bool IsReadOnly
        {
            get { return true; }
        }

        #endregion

        #region IList

        public int IndexOf(TOutputType item)
        {
            throw new NotImplementedException();
        }

        public void Insert(int index, TOutputType item)
        {
            throw new NotImplementedException();
        }

        public void RemoveAt(int index)
        {
            throw new NotImplementedException();
        }

        public TOutputType this[int index]
        {
            get
            {
                if (index < 0)
                    throw new IndexOutOfRangeException();

                int begin = 0;
                int end = 0;

                foreach (DataSource<TOutputType> p in parts)
                {
                    end = end + p.Count;
                    if (index < end)
                        return p[index - begin];

                    begin = end;
                }

                throw new IndexOutOfRangeException();
            }

            set
            {
                // OK, not implemented in ReaderBase
                throw new NotImplementedException();
            }
        }

        #endregion
    }
}
