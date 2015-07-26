using System;
using Bio;
using Bio.IO.BAM;
using Bio.IO.SAM;
using System.Collections.Generic;
using System.Linq;

namespace Bio.IO.PacBio
{
    public class PacBioBamReader
    {
        public static IEnumerable<PacBioCCSRead> ParseReads(string fileName) {
            var bp = new BAMParser ();
            var seqs = bp.Parse (fileName);
            foreach (var s in seqs.AlignedSequences) {
                yield return new PacBioCCSRead (s as SAMAlignedSequence);
            }
        }
    }
}

