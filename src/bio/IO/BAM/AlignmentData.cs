using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Runtime.InteropServices;
using Bio.IO.SAM;

namespace Bio.IO.BAM
{
    /// <summary>
    /// A BAM header section, useful for de-serializing the SAM alignment header 
    /// data with one pointer cast instead of unpacking bits.
    /// 
    /// ACHTUNG!!!:  The data are stored in little endian format, so the items 
    /// marked in the spec as the upper end of the << shifts are actually at 
    /// the bottom of the int.
    /// </summary>
    [StructLayout(LayoutKind.Sequential)]
    struct AlignmentData
    {
        //int 1 section
        public int refSeqIndex;
        //int 2 section
        public int alignedSeqPosition;
        //int 3 sectionspec:  bin<<16|MAPQ<<8|l read name;
        public byte ReadNameLength;
        public byte MAPQ;
        public ushort bin;
        //int 4 section spec: FLAG<<16|n cigar op;
        public ushort cigarLen;
        public SAMFlags flagValue;//size ushort
        //int 5 section
        public int readLen;
        //int 6 section
        public int mateRefSeqIndex;
        //int 7 section
        public int mateSeqPosition;
        //int 8 section
        public int templateLength;
    }
}
