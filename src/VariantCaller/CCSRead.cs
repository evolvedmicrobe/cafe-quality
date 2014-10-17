using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using Bio;

namespace VariantCaller
{
    public class CCSRead
    {
        public readonly int ZMW;
        public List<CCSSubRead> SubReads;
        public byte[] Seq;
        
        /// <summary>
        /// The DNA this read is derived from.
        /// </summary>
        public Reference AssignedRefrerence;

        public CCSRead(Sequence read)
        {
            ZMW = Convert.ToInt32(read.ID.Split('/')[1]);
            Seq = read.ToArray();
        } 
    }
}
