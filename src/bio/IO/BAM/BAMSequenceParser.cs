//#define WANT_OLD_VERSION
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Bio.IO.BAM;
using Bio;
using System.IO;
using Bio.IO.SAM;
using Bio.Util;
using System.Globalization;
using System.Diagnostics;
using Bio.IO;

namespace Bio.IO.BAM
{
    /// <summary>
    /// An implementation of a BAM parser that only returns sequences, not full alignments
    /// </summary>
    public class BAMSequenceParser : BAMParser
    {
        private string _fileName;
        public string ChromosomeToGet;

        public BAMSequenceParser(string Filename)
        {
            this.Open(Filename);
        }
        public void Open(string dataSource)
        {
            this._fileName = dataSource;
        }

#if WANT_OLD_VERSION
         /// <summary>
        /// Returns an aligned sequence by parses the BAM file.
        /// </summary>
        private SAMAlignedSequence GetAlignedSequence(int start, int end)
        {
            byte[] array = new byte[4];

            ReadUnCompressedData(array, 0, 4);
            int blockLen = Helper.GetInt32(array, 0);
            byte[] alignmentBlock = new byte[blockLen];
            ReadUnCompressedData(alignmentBlock, 0, blockLen);
            SAMAlignedSequence alignedSeq = new SAMAlignedSequence();
            int value;
            UInt32 UnsignedValue;
            // 0-4 bytes
            int refSeqIndex = Helper.GetInt32(alignmentBlock, 0);

            if (refSeqIndex == -1)
                alignedSeq.RName = "*";
            else
                alignedSeq.RName = refSeqNames[refSeqIndex];
            

            // 4-8 bytes
            alignedSeq.Pos = Helper.GetInt32(alignmentBlock, 4) + 1;

            // if there is no overlap no need to parse further.
            //     BAMPos > closedEnd
            // => (alignedSeq.Pos - 1) > end -1
            if (alignedSeq.Pos > end)
            {
                return null;
            }

            // 8 - 12 bytes "bin<<16|mapQual<<8|read_name_len"
            UnsignedValue = Helper.GetUInt32(alignmentBlock, 8);

            // 10 -12 bytes
            alignedSeq.Bin = (int)(UnsignedValue & 0xFFFF0000) >> 16;
            // 9th bytes
            alignedSeq.MapQ = (int)(UnsignedValue & 0x0000FF00) >> 8;
            // 8th bytes
            int queryNameLen = (int)(UnsignedValue & 0x000000FF);

            // 12 - 16 bytes
            UnsignedValue = Helper.GetUInt32(alignmentBlock, 12);
            // 14-16 bytes
            int flagValue = (int)(UnsignedValue & 0xFFFF0000) >> 16;
            alignedSeq.Flag = (SAMFlags)flagValue;
            // 12-14 bytes
            int cigarLen = (int)(UnsignedValue & 0x0000FFFF);

            // 16-20 bytes
            int readLen = Helper.GetInt32(alignmentBlock, 16);

            // 20-24 bytes
            int mateRefSeqIndex = Helper.GetInt32(alignmentBlock, 20);
            if (mateRefSeqIndex != -1)
            {
                alignedSeq.MRNM = refSeqNames[mateRefSeqIndex];
            }
            else
            {
                alignedSeq.MRNM = "*";
            }

            // 24-28 bytes
            alignedSeq.MPos = Helper.GetInt32(alignmentBlock, 24) + 1;

            // 28-32 bytes
            alignedSeq.ISize = Helper.GetInt32(alignmentBlock, 28);

            // 32-(32+readLen) bytes
            alignedSeq.QName = System.Text.ASCIIEncoding.ASCII.GetString(alignmentBlock, 32, queryNameLen - 1);
            StringBuilder strbuilder = new StringBuilder();
            int startIndex = 32 + queryNameLen;

            for (int i = startIndex; i < (startIndex + cigarLen * 4); i += 4)
            {
                // Get the CIGAR operation length stored in first 28 bits.
                UInt32 cigarValue = Helper.GetUInt32(alignmentBlock, i);
                strbuilder.Append(((cigarValue & 0xFFFFFFF0) >> 4).ToString(CultureInfo.InvariantCulture));

                // Get the CIGAR operation stored in last 4 bits.
                value = (int)cigarValue & 0x0000000F;

                // MIDNSHP=>0123456
                switch (value)
                {
                    case 0:
                        strbuilder.Append("M");
                        break;
                    case 1:
                        strbuilder.Append("I");
                        break;
                    case 2:
                        strbuilder.Append("D");
                        break;
                    case 3:
                        strbuilder.Append("N");
                        break;
                    case 4:
                        strbuilder.Append("S");
                        break;
                    case 5:
                        strbuilder.Append("H");
                        break;
                    case 6:
                        strbuilder.Append("P");
                        break;
                    case 7:
                        strbuilder.Append("=");
                        break;
                    case 8:
                        strbuilder.Append("X");
                        break;
                    default:
                        throw new FileFormatException(Properties.Resource.BAM_InvalidCIGAR);
                }
            }

            string cigar = strbuilder.ToString();
            if (string.IsNullOrWhiteSpace(cigar))
            {
                alignedSeq.CIGAR = "*";
            }
            else
            {
                alignedSeq.CIGAR = cigar;
            }

            // if there is no overlap no need to parse further.
            // ZeroBasedRefEnd < start
            // => (alignedSeq.RefEndPos -1) < start
            if (alignedSeq.RefEndPos - 1 < start && alignedSeq.RName!=Properties.Resource.SAM_NO_REFERENCE_DEFINED_INDICATOR)
            {
                return null;
            }

            startIndex += cigarLen * 4;
            strbuilder = new StringBuilder();
            int index = startIndex;
            for (; index < (startIndex + (readLen + 1) / 2) - 1; index++)
            {
                // Get first 4 bit value
                value = (alignmentBlock[index] & 0xF0) >> 4;
                strbuilder.Append(GetSeqChar(value));
                // Get last 4 bit value
                value = alignmentBlock[index] & 0x0F;
                strbuilder.Append(GetSeqChar(value));
            }

            value = (alignmentBlock[index] & 0xF0) >> 4;
            strbuilder.Append(GetSeqChar(value));
            if (readLen % 2 == 0)
            {
                value = alignmentBlock[index] & 0x0F;
                strbuilder.Append(GetSeqChar(value));
            }

            startIndex = index + 1;
            string strSequence = strbuilder.ToString();
            byte[] qualValues = new byte[readLen];
            string strQualValues = "*";

            if (alignmentBlock[startIndex] != 0xFF)
            {
                for (int i = startIndex; i < (startIndex + readLen); i++)
                {
                    qualValues[i - startIndex] = (byte)(alignmentBlock[i] + 33);
                }

                strQualValues = System.Text.ASCIIEncoding.ASCII.GetString(qualValues);
            }

            SAMParser.ParseQualityNSequence(alignedSeq, Alphabet, strSequence, strQualValues);

            startIndex += readLen;
           
            if (alignmentBlock.Length > startIndex + 4 && alignmentBlock[startIndex] != 0x0 && alignmentBlock[startIndex + 1] != 0x0)
            {

                for (index = startIndex; index < alignmentBlock.Length; )
                {
                    SAMOptionalField optionalField = new SAMOptionalField();
                    optionalField.Tag = System.Text.ASCIIEncoding.ASCII.GetString(alignmentBlock, index, 2);
                    index += 2;
                    char vType = (char)alignmentBlock[index++];
                    string valueType = vType.ToString();

                    // SAM format supports [AifZH] for value type.
                    // In BAM, an integer may be stored as a signed 8-bit integer (c), unsigned 8-bit integer (C), signed short (s), unsigned
                    // short (S), signed 32-bit (i) or unsigned 32-bit integer (I), depending on the signed magnitude of the integer. However,
                    // in SAM, all types of integers are presented as type ʻiʼ. 
                    string message = Helper.IsValidPatternValue("VType", valueType, BAMOptionalFieldRegex);
                    if (!string.IsNullOrEmpty(message))
                    {
                        throw new FormatException(message);
                    }


                    optionalField.Value = GetOptionalValue(vType, alignmentBlock, ref index).ToString();

                    // Convert to SAM format.
                    if ("cCsSI".IndexOf(vType) >= 0)
                    {
                        valueType = "i";
                    }

                    optionalField.VType = valueType;

                    alignedSeq.OptionalFields.Add(optionalField);
                }
            }

            return alignedSeq;
        }
#endif
#if !WANT_OLD_VERSION

		protected CompactSAMSequence GetAlignedSequence()
        {
            byte[] array = new byte[4];
            ReadUnCompressedData(array, 0, 4);
            int blockLen = Helper.GetInt32(array, 0);
            byte[] alignmentBlock = new byte[blockLen];
            ReadUnCompressedData(alignmentBlock, 0, blockLen);
            int value;
            UInt32 UnsignedValue;
            // 0-4 bytes
            int refSeqIndex = Helper.GetInt32(alignmentBlock, 0);
            string RName;
            if (refSeqIndex == -1)
                RName = "*";
            else
                RName = refSeqNames[refSeqIndex];

            // 4-8 bytes
            int Pos = Helper.GetInt32(alignmentBlock, 4) + 1;

            // 8 - 12 bytes "bin<<16|mapQual<<8|read_name_len"
            UnsignedValue = Helper.GetUInt32(alignmentBlock, 8);
            int queryNameLen = (int)(UnsignedValue & 0x000000FF);

            // 12 - 16 bytes
            UnsignedValue = Helper.GetUInt32(alignmentBlock, 12);
            int flagValue = (int)(UnsignedValue & 0xFFFF0000) >> 16;
            int cigarLen = (int)(UnsignedValue & 0x0000FFFF);

            //// 16-20 bytes
            int readLen = Helper.GetInt32(alignmentBlock, 16);

            // 32-(32+readLen) bytes
            string name = System.Text.ASCIIEncoding.ASCII.GetString(alignmentBlock, 32, queryNameLen - 1);
            StringBuilder strbuilder = new StringBuilder();
            int startIndex = 32 + queryNameLen;
            for (int i = startIndex; i < (startIndex + cigarLen * 4); i += 4)
            {
                // Get the CIGAR operation length stored in first 28 bits.
                UInt32 cigarValue = Helper.GetUInt32(alignmentBlock, i);
                strbuilder.Append(((cigarValue & 0xFFFFFFF0) >> 4).ToString(CultureInfo.InvariantCulture));

                // Get the CIGAR operation stored in last 4 bits.
                value = (int)cigarValue & 0x0000000F;

                // MIDNSHP=>0123456
                switch (value)
                {
                    case 0:
                        strbuilder.Append("M");
                        break;
                    case 1:
                        strbuilder.Append("I");
                        break;
                    case 2:
                        strbuilder.Append("D");
                        break;
                    case 3:
                        strbuilder.Append("N");
                        break;
                    case 4:
                        strbuilder.Append("S");
                        break;
                    case 5:
                        strbuilder.Append("H");
                        break;
                    case 6:
                        strbuilder.Append("P");
                        break;
                    case 7:
                        strbuilder.Append("=");
                        break;
                    case 8:
                        strbuilder.Append("X");
                        break;
                    default:
                        throw new FileFormatException(Properties.Resource.BAM_InvalidCIGAR);
                }
            }

            string cigar = strbuilder.ToString();
            if (string.IsNullOrWhiteSpace(cigar))
            {
                cigar = "*";
            }
            startIndex += cigarLen * 4;
            //strbuilder = new StringBuilder();
            byte[] seqData = new byte[readLen];
            int seqDataIndex = 0;
            int index = startIndex;
            for (; index < (startIndex + (readLen + 1) / 2) - 1; index++)
            {
                // Get first 4 bit value
                value = (alignmentBlock[index] & 0xF0) >> 4;
                //strbuilder.Append(GetSeqChar(value));
                seqData[seqDataIndex++] = GetSeqCharAsByte(value);
                // Get last 4 bit value
                value = alignmentBlock[index] & 0x0F;
                //strbuilder.Append(GetSeqChar(value));
                seqData[seqDataIndex++]=GetSeqCharAsByte(value);
            }
            value = (alignmentBlock[index] & 0xF0) >> 4;
            //strbuilder.Append(GetSeqChar(value));
            seqData[seqDataIndex++] = GetSeqCharAsByte(value);
            if (readLen % 2 == 0)
            {
                value = alignmentBlock[index] & 0x0F;
                //strbuilder.Append(GetSeqChar(value));
                seqData[seqDataIndex++] = GetSeqCharAsByte(value);
            }
            startIndex = index + 1;
           // string strSequence = strbuilder.ToString();
            //Insert qual value catch here?  ADDING NEW QUALITY SCORE FINDER!!!
            byte[] qualValues = new byte[readLen];
            string strQualValues = "*";
            if (alignmentBlock[startIndex] != 0xFF)
            {
                for (int i = startIndex; i < (startIndex + readLen); i++)
                {
                    qualValues[i - startIndex] = (byte)(alignmentBlock[i] + 33);
                }

                strQualValues = System.Text.ASCIIEncoding.ASCII.GetString(qualValues);
            }
            //END NEW EDITION!

            //var syms = Encoding.UTF8.GetBytes(strSequence);
            var alpha = Alphabets.AutoDetectAlphabet(seqData, 0, seqData.Length, null);
            //Sequence toReturn = new Sequence(alpha, syms);
            //TODO: Possibly a bit unsafe here
            var toReturn = new CompactSAMSequence(alpha, FastQFormatType.GATK_Recalibrated, seqData, qualValues, false);
            toReturn.ID = name;
            toReturn.Pos = Pos;
            toReturn.CIGAR = cigar;
            toReturn.RName = RName;
			toReturn.SAMFlags = (SAMFlags)flagValue;
            return toReturn;

        }
      
#endif

        public IEnumerable<CompactSAMSequence> Parse()
        {
            if (string.IsNullOrWhiteSpace(_fileName))
            {
                throw new ArgumentNullException("fileName");
            }
            using (readStream = new FileStream(_fileName, FileMode.Open, FileAccess.Read, FileShare.Read))
            {
                Stream reader = readStream;
                if (reader == null || reader.Length == 0)
                {
                    throw new FileFormatException(Properties.Resource.BAM_InvalidBAMFile);
                }
                if (!String.IsNullOrEmpty(ChromosomeToGet))
                {
                    foreach (var s in
                    ParseRangeAsEnumerableSequences(_fileName, ChromosomeToGet ))
                    {
						if (s != null) {
							yield return s;
						}
                        ////TODO: Super inefficient right now, am parsing the sequence multiple times,
                        ////fix this.
                        //var s2 = s.ToArray ();
                        //var alpha = Alphabets.AutoDetectAlphabet(s2, 0, s2.Length, null);


                        //var strippedOfInfo = new Sequence(alpha, s2);
                        //yield return strippedOfInfo;
                    }
                }
                else
                {
                    readStream = reader;
                    ValidateReader();
                    SAMAlignmentHeader header = GetHeader();
                    SequenceAlignmentMap sequenceAlignmentMap = null;
                    if (sequenceAlignmentMap == null)
                    {
                        sequenceAlignmentMap = new SequenceAlignmentMap(header);
                    }

                    while (!IsEOF())
                    {

#if WANT_OLD_VERSION
                    SAMAlignedSequence alignedSeq = GetAlignedSequence(0, int.MaxValue);
#else
						var alignedSeq = GetAlignedSequence();
#endif
                        if (alignedSeq != null)
                        {

#if WANT_OLD_VERSION
                            //make a new Sequence
                            ISequence strippedOfInfo = null;
                            try
                            {
                            var syms=alignedSeq.QuerySequence.ToArray();
                            var alpha = Alphabets.AutoDetectAlphabet(syms, 0, syms.Length, null);
                            strippedOfInfo = new Sequence(alpha, alignedSeq.QuerySequence.ToArray());

                                strippedOfInfo = alignedSeq;

                            }
                            catch (ArgumentOutOfRangeException exception)
                            {
                                Debug.Write("Could not convert sequence: " + exception.Message);
                            }
                            if (strippedOfInfo != null)
                                yield return strippedOfInfo;
#else
							yield return alignedSeq;
							#endif
                        }
                        alignedSeq = null;
                    }
                }
            }
        }
        public IEnumerable<ISequence> Parse(System.IO.StreamReader reader)
        {
            throw new NotImplementedException();
        }
        public SAMAlignmentHeader GetFileHeader() 
        {
            SAMAlignmentHeader header;
            using (FileStream bamStream = new FileStream(_fileName, FileMode.Open, FileAccess.Read, FileShare.Read))
            {
                string bamIndexFileName = getBAMIndexFileName(_fileName);
                using (BAMIndexFile bamIndexFile = new BAMIndexFile(bamIndexFileName, FileMode.Open, FileAccess.Read))
                {
                    readStream = bamStream;
                    if (readStream == null || readStream.Length == 0)
                    {
                        throw new FileFormatException(Properties.Resource.BAM_InvalidBAMFile);
                    }
                    ValidateReader();
                    header = GetHeader();
                }
            }
            return header;
        }

		public IEnumerable<CompactSAMSequence> ParseRangeAsEnumerableSequences(string fileName, string refSeqName, int start = 0, int end = Int32.MaxValue)
        {
            if (refSeqName == null)
            {
                throw new ArgumentNullException("refSeqName");
            }
            using (FileStream bamStream = new FileStream(fileName, FileMode.Open, FileAccess.Read, FileShare.Read))
            {
                string bamIndexFileName = getBAMIndexFileName(fileName);
                using (BAMIndexFile bamIndexFile = new BAMIndexFile(bamIndexFileName, FileMode.Open, FileAccess.Read))
                {
                    readStream = bamStream;
                    if (readStream == null || readStream.Length == 0)
                    {
                        throw new FileFormatException(Properties.Resource.BAM_InvalidBAMFile);
                    }
                    ValidateReader();
                    SAMAlignmentHeader header = GetHeader();
                    // verify whether there is any reads related to chromosome.
                    int refSeqIndex = refSeqNames.IndexOf(refSeqName);
                    if (refSeqIndex < 0)
                    {
                        string message = string.Format(CultureInfo.InvariantCulture, Properties.Resource.BAM_RefSeqNotFound, refSeqName);
                        throw new ArgumentException(message, "refSeqName");
                    }
                    BAMIndex bamIndexInfo = bamIndexFile.Read();
                    BAMReferenceIndexes refIndex = bamIndexInfo.RefIndexes[refSeqIndex];
                    IList<Chunk> chunks = GetChunks(refIndex, start, end);
                    foreach (var s in EnumerateAlignedSequences(chunks))
                    {
                        if (s != null &&  (s.RName == "*" || (s.Pos >= (start - 1) && s.RefEndPos < end) ))
                        {
                            yield return s;
                        }
                    }
                    readStream = null;
                }
            }
        }
		private IEnumerable<CompactSAMSequence> EnumerateAlignedSequences(IList<Chunk> chunks)
        {
            foreach (Chunk chunk in chunks)
            {
                readStream.Seek((long)chunk.ChunkStart.CompressedBlockOffset, SeekOrigin.Begin);
                GetNextBlock();
                if (deCompressedStream != null)
                {
                    deCompressedStream.Seek(chunk.ChunkStart.UncompressedBlockOffset, SeekOrigin.Begin);
                    // read until eof or end of the chunck is reached.
                    while (!IsEOF() && (currentCompressedBlockStartPos < (long)chunk.ChunkEnd.CompressedBlockOffset || deCompressedStream.Position < chunk.ChunkEnd.UncompressedBlockOffset))
                    {
                        var alignedSeq = GetAlignedSequence();
                        if (alignedSeq != null)
                        {
                            yield return alignedSeq;
                        }
                    }
                }
            }
        }
        public void Close()
        {

        }
        public IAlphabet Alphabet
        {
            get
            {
                return base.Alphabet;
            }
            set
            {
                throw new NotSupportedException(Properties.Resource.BAMParserAlphabetCantBeSet);
            }
        }

        public string Name
        {
            get { return "BAM Sequence Parser"; }
        }

        public string Description
        {
            get { throw new NotImplementedException(); }
        }

        public string SupportedFileTypes
        {
            get { return ".bam"; }
        }

        public void Dispose()
        {
        }
    }
}
