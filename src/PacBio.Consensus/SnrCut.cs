using System;
using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Diagnostics;
using System.Linq;
using ConsensusCore;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics.LinearAlgebra;

using PacBio.Utils;

using Vector=MathNet.Numerics.LinearAlgebra.Double.DenseVector;


namespace PacBio.Consensus
{
    // (more general implementation is coming in PacBio.Data)
    internal struct Range
    {
        public readonly float Begin;
        public readonly float End;

        public Range(float begin, float end) 
        {
            Begin = begin;
            End = end;
        }

        public static readonly Range Universe = new Range(float.MinValue, float.MaxValue);

        public bool Contains(float pt)
        {
            return (Begin <= pt && pt < End);
        }

        public static Range Parse(string s)
        {
            var elts = s.Split(new char[] {'-'});
            return new Range(float.Parse(elts[0]), float.Parse(elts[1]));
        }
    }
    public class SnrCut
    {
     
            // in TGAC order
            private Range[] Ranges;
            private SnrCut() {}

            public static readonly SnrCut PassAll = new SnrCut { Ranges =  4.Fill(Range.Universe) };

            private static Tuple<char, Range> MakeRange(string clause)
            {   
                char channel = clause[0];
                var range = Range.Parse(clause.Substring(1));
                return Tuple.Create(channel, range);
            }

            // format expected: A3-5,T6-8
            // omitted channels have no filtering applied
            public static SnrCut Parse(string s)
            {
                var clauses = s.Split(new char[] {','}, StringSplitOptions.RemoveEmptyEntries);
                var ranges = 4.Fill(Range.Universe);
                clauses.Select(MakeRange).ForEach(t =>
                    {
                        switch(t.Item1) {
                        case 'T': ranges[0] = t.Item2; break;
                        case 'G': ranges[1] = t.Item2; break;
                        case 'A': ranges[2] = t.Item2; break;
                        case 'C': ranges[3] = t.Item2; break;
                        }
                    });
                return new SnrCut { Ranges = ranges };
            }

            public bool Contains(float[] snrs)
            {
                return Ranges.Zip(snrs, (r, snr) => r.Contains(snr)).All();
            }
    }

}

