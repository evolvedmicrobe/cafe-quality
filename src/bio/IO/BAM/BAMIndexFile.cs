﻿using System;
using System.IO;
using Bio.Util;
using System.Collections;

namespace Bio.IO.BAM
{
    /// <summary>
    /// Class to read or write BAMIndex data from a file or a stream.
    /// </summary>
    public class BAMIndexFile : IDisposable
    {
        #region Constants and Static Methods
        /// <summary>
        /// The highest number of bins allowed, meta-data can be stored in the chunks position for a bin
        /// this large
        /// </summary>
        internal const int MAX_BINS = 37450;   // =(8^6-1)/7+1

        /// <summary>
        /// The number of 16kb (2^14) bins in the indexing scheme
        /// </summary>
        internal const int MAX_LINERINDEX_ARRAY_SIZE=MAX_BINS+1-4681;

        /// <summary>
        /// Not all sequences can get all possible bins, so this returns the largest sequence length possible
        /// </summary>
        /// <param name="sequenceLength"></param>
        /// <returns></returns>
        internal static int LargestBinPossibleForSequenceLength(int sequenceLength)
        {
            return 4681 + (sequenceLength >> 14);
        }

        #endregion 
        #region Private Fields
        // holds index stream
        private Stream sourceStream;
        #endregion

        #region Properties
        /// <summary>
        /// Gets the underlying stream.
        /// </summary>
        public Stream Source
        {
            get
            {
                return sourceStream;
            }
        }
        #endregion

        #region Constructors and Destructors
        /// <summary>
        /// Creates new instance of BAMIndexFile class with specified filename.
        /// </summary>
        /// <param name="filename">Index filename to use while reading or writing BAMIndex data.</param>
        /// <param name="mode">File mode to use while creating or opening specified file.</param>
        /// <param name="access">File access to use while creating or opening specified file.</param>
        public BAMIndexFile(string filename, FileMode mode, FileAccess access)
        {
            if (string.IsNullOrWhiteSpace(filename))
            {
                throw new ArgumentNullException("filename");
            }

            sourceStream = new FileStream(filename, mode, access);
        }

        /// <summary>
        /// Creates new instance of BAMIndexFile clas with specified stream.
        /// </summary>
        /// <param name="stream">Stream to use while reading or writing BAMIndex data.</param>
        public BAMIndexFile(Stream stream)
        {
            if (stream == null)
            {
                throw new ArgumentNullException("stream");
            }

            sourceStream = stream;
        }
        #endregion

        #region Public Methods
        /// <summary>
        /// Writes specified BAMIndex data.
        /// </summary>
        /// <param name="bamIndex">BAMIndex instance to write.</param>
        public void Write(BAMIndex bamIndex)
        {
            if (bamIndex == null)
            {
                throw new ArgumentNullException("bamIndex");
            }
            if (sourceStream == null)
            {
                throw new InvalidOperationException(Properties.Resource.BAM_CantUseBAMIndexStreamDisposed);
            }
            byte[] arrays = new byte[20];

            byte[] magic = new byte[] { 66, 65, 73, 1 };
            Write(magic, 0, 4);

            arrays = Helper.GetLittleEndianByteArray(bamIndex.RefIndexes.Count);
            Write(arrays, 0, 4);

            for (Int32 refindex = 0; refindex < bamIndex.RefIndexes.Count; refindex++)
            {
                
                BAMReferenceIndexes bamindices = bamIndex.RefIndexes[refindex];
                int binCount = bamindices.Bins.Count;
                bool addingMetaData = bamindices.HasMetaData && BitConverter.IsLittleEndian;
                if (addingMetaData)
                {
                    binCount++;
                }                
                arrays = Helper.GetLittleEndianByteArray(binCount);
                Write(arrays, 0, 4);
                //Write each bin
                for (Int32 binIndex = 0; binIndex < bamindices.Bins.Count; binIndex++)
                {
                    Bin bin = bamindices.Bins[binIndex];
                    arrays = Helper.GetLittleEndianByteArray(bin.BinNumber);
                    Write(arrays, 0, 4);
                    int chunkCount = bin.Chunks.Count;
                   
                    arrays = Helper.GetLittleEndianByteArray(chunkCount);
                    Write(arrays, 0, 4);
                    for (Int32 chunkIndex = 0; chunkIndex < bin.Chunks.Count; chunkIndex++)
                    {
                        Chunk chunk = bin.Chunks[chunkIndex];
                        arrays = GetBAMOffsetArray(chunk.ChunkStart);
                        Write(arrays, 0, 8);
                        arrays = GetBAMOffsetArray(chunk.ChunkEnd);
                        Write(arrays, 0, 8);
                    }
                    
                }
                //Add Meta Data - this varies by implementation, .NET Bio will do start and
                //end of reads found in file and then mapped/unmapped
                //TODO: Assumes little endian, only adds if so
                if (addingMetaData)
                {
                    //Dummy bin to indicate meta-data
                    arrays = Helper.GetLittleEndianByteArray(BAMIndexFile.MAX_BINS);
                    Write(arrays, 0, 4);
                    //2 chunks worth of meta data
                    //first the file offsets
                    arrays = Helper.GetLittleEndianByteArray((int)2);
                    Write(arrays, 0, 4);
                    arrays = GetBAMOffsetArray(bamindices.FirstOffSetSeen);
                    Write(arrays, 0, 8);
                    arrays = GetBAMOffsetArray(bamindices.LastOffSetSeen);
                    Write(arrays, 0, 8);
                    arrays = BitConverter.GetBytes(bamindices.MappedReadsCount);
                    Write(arrays, 0, 8);
                    arrays = BitConverter.GetBytes(bamindices.UnMappedReadsCount);
                    Write(arrays, 0, 8);
                }
                arrays = Helper.GetLittleEndianByteArray(bamindices.LinearIndex.Count);
                Write(arrays, 0, 4);
                for (Int32 offsetIndex = 0; offsetIndex < bamindices.LinearIndex.Count; offsetIndex++)
                {
                    FileOffset value = bamindices.LinearIndex[offsetIndex];
                    arrays = GetBAMOffsetArray(value);
                    Write(arrays, 0, 8);
                }
                sourceStream.Flush();
            }
            sourceStream.Flush();
        }

        /// <summary>
        /// Returns BAMIndex instance by parsing BAM index source.
        /// </summary>
        public BAMIndex Read()
        {
            if (sourceStream == null)
            {
                throw new InvalidOperationException(Properties.Resource.BAM_CantUseBAMIndexStreamDisposed);
            }

            BAMIndex bamIndex = new BAMIndex();
            byte[] arrays = new byte[20];

            Read(arrays, 0, 4);

            if (arrays[0] != 66 || arrays[1] != 65 || arrays[2] != 73 || arrays[3] != 1)
            {
                throw new FormatException(Properties.Resource.BAM_InvalidIndexFile);
            }
            Read(arrays, 0, 4);
            int n_ref = Helper.GetInt32(arrays, 0);
            for (Int32 refindex = 0; refindex < n_ref; refindex++)
            {
                BAMReferenceIndexes bamindices = new BAMReferenceIndexes();
                bamIndex.RefIndexes.Add(bamindices);
                Read(arrays, 0, 4);
                int n_bin = Helper.GetInt32(arrays, 0);
                for (Int32 binIndex = 0; binIndex < n_bin; binIndex++)
                {
                    Bin bin = new Bin();
                    Read(arrays, 0, 4);
                    bin.BinNumber = Helper.GetUInt32(arrays, 0);
                    Read(arrays, 0, 4);
                    int n_chunk = Helper.GetInt32(arrays, 0);
                    if (bin.BinNumber == MAX_BINS)//some groups use this to place meta-data, such as the picard toolkit and now SAMTools
                    {
                        //Meta data was later added in to the SAMTools specification
                        for (Int32 chunkIndex = 0; chunkIndex < n_chunk; chunkIndex++)
                        {
                            bamindices.HasMetaData = true;
                            Read(arrays, 0, 8);
                            bamindices.MappedReadsCount = Helper.GetUInt64(arrays, 0);
                            Read(arrays, 0, 8);
                            bamindices.UnMappedReadsCount = Helper.GetUInt64(arrays, 0);
                        }

                    }
                    else if (bin.BinNumber > MAX_BINS)
                    {
                        throw new FileFormatException("BAM Index is incorrectly formatted.  Bin number specified is higher than the maximum allowed.");
                    }
                    else
                    {
                         bamindices.Bins.Add(bin);
                        for (Int32 chunkIndex = 0; chunkIndex < n_chunk; chunkIndex++)
                        {
                            Chunk chunk = new Chunk();
                            bin.Chunks.Add(chunk);
                            Read(arrays, 0, 8);
                            chunk.ChunkStart = GetBAMOffset(arrays, 0);
                            Read(arrays, 0, 8);
                            chunk.ChunkEnd = GetBAMOffset(arrays, 0);
                        }
                    }
                }
                //Get number of linear bins
                Read(arrays, 0, 4);
                int n_intv = Helper.GetInt32(arrays, 0);

                for (Int32 offsetIndex = 0; offsetIndex < n_intv; offsetIndex++)
                {
                    FileOffset value;
                    Read(arrays, 0, 8);
                    value = GetBAMOffset(arrays, 0);
                    bamindices.LinearIndex.Add(value);
                }
            }
            

            return bamIndex;
        }
        #endregion

        #region IDisposable Members
        /// <summary>
        /// Disposes resources held by this object.
        /// </summary>
        public void Dispose()
        {
            Dispose(true);
            GC.SuppressFinalize(this);
        }

        /// <summary>
        /// Dispose the underlying stream.
        /// </summary>
        /// <param name="disposing">If disposing equals true, Requests that the system not call the finalizer for this instance.</param>
        protected virtual void Dispose(bool disposing)
        {
            if (sourceStream != null)
            {
                sourceStream.Close();
                sourceStream.Dispose();
                sourceStream = null;
            }
        }

        #endregion

        #region Private Methods
        // Converts bytes array to FileOffset object.
        private static FileOffset GetBAMOffset(byte[] bytes, int startIndex)
        {
            
            UInt64 value = bytes[startIndex + 7];
            value = (value << 8) + bytes[startIndex + 6];
            value = (value << 8) + bytes[startIndex + 5];
            value = (value << 8) + bytes[startIndex + 4];
            value = (value << 8) + bytes[startIndex + 3];
            value = (value << 8) + bytes[startIndex + 2];
            UInt16 uvalue = bytes[startIndex + 1];
            uvalue = (UInt16)((UInt16)(uvalue << 8) + (UInt16)bytes[startIndex]);
            FileOffset offset = new FileOffset(value,uvalue);
            return offset;
        }

        // Converts FileOffset object to byte array.
        private static byte[] GetBAMOffsetArray(FileOffset offset)
        {
            byte[] bytes = new byte[8];

            bytes[0] = (byte)(offset.UncompressedBlockOffset & 0x00FF);
            bytes[1] = (byte)((offset.UncompressedBlockOffset & 0xFF00) >> 8);

            bytes[2] = (byte)(offset.CompressedBlockOffset & 0x0000000000FF);
            bytes[3] = (byte)((offset.CompressedBlockOffset & 0x00000000FF00) >> 8);
            bytes[4] = (byte)((offset.CompressedBlockOffset & 0x000000FF0000) >> 16);
            bytes[5] = (byte)((offset.CompressedBlockOffset & 0x0000FF000000) >> 24);
            bytes[6] = (byte)((offset.CompressedBlockOffset & 0x00FF00000000) >> 32);
            bytes[7] = (byte)((offset.CompressedBlockOffset & 0xFF0000000000) >> 40);

            return bytes;
        }

        // Writes byte array to underlying stream of this instance.
        private void Write(byte[] array, int offset, int count)
        {
            sourceStream.Write(array, offset, count);
        }

        // reads specified number of bytes from the underlying stream to specified array starting from specified offset.
        private void Read(byte[] array, int offset, int count)
        {
            if (IsEOF() || sourceStream.Read(array, offset, count) != count)
            {
                throw new FileFormatException(Properties.Resource.BAM_InvalidIndexFile+"\nCould not read the correct number of bytes from file");
            }
        }

        // Gets a boolean which indicates whether underlying stream reached EOF or not.
        private bool IsEOF()
        {
            return sourceStream.Position == sourceStream.Length;
        }
        #endregion

    }
}
