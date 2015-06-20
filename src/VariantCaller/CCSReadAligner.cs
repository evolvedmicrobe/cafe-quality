using System;
using System.Collections.Generic;
using VariantCaller;
using System.Linq;
using PacBio.Utils;

namespace VariantCaller
{
    public class CCSReadAligner
    {
        public static Tuple<string, List<Variant>>  CallVariants(CCSRead seq)
        {   
            if(seq.AssignedReference == null) {return null;}
            var alns = seq.AssignedReference.AlignSequence(seq.Seq).ToArray();
            if (alns.Length ==0 ) {
                return null;
            }
            var best_score = alns.Max (z => z.Score);
            var best = alns.First (z => z.Score == best_score);  
            var variants = VariantCaller.CallVariants (best, seq.AssignedReference.RefSeq, seq.Seq);
            return new Tuple<string, List<Variant>>(best.ToString(), variants);            
        }
    }
}

