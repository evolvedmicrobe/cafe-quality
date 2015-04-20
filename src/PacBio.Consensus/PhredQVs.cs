using System;
using System.Collections.Generic;

namespace PacBio.Consensus
{
    public static class QVs
    {
        public static double CombineErrorProbability(IEnumerable<double> d)
        {
            var noErrorProb = 1.0;

            foreach (var errorProb in d) 
            {
                noErrorProb *= (1.0 - errorProb);
            }

            return 1.0 - noErrorProb;
        }

        public static double CombineErrorProbability(params double[] d)
        {
            var noErrorProb = 1.0;

            foreach (var errorProb in d) 
            {
                noErrorProb *= (1.0 - errorProb);
            }

            return 1.0 - noErrorProb;
        }

        public static double PhredProb(byte qv)
        {
            return Math.Pow(10, qv / -10.0);
        }

        public static double PhredProb(double qv)
        {
            return Math.Pow(10, qv / -10.0);
        }

        public static byte PhredQV(double pErr)
        {
            pErr = Math.Min(0.99, pErr);
            return (byte)Math.Max(0,  Math.Round(-10 * Math.Log10(pErr)));
        }

        public static byte ProbToQV(double pErr)
        {
            return PhredQV(pErr);
        }

        public static double QVToProb(byte qv)
        {
            return PhredProb(qv);
        }
    }
}

