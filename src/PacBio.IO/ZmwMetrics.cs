using System;
using System.Linq;
using PacBio.HDF;
using PacBio.Utils;

namespace PacBio.IO
{
    /// <summary>
    /// A static class to compute ZMW metrics that are derived from ZmwBases data
    /// </summary>
    public static class ZmwMetricsBasesFunc
    {
        private static int nChan = 4;

        /// <summary>
        /// Return the fraction of calls made by channel
        /// </summary>
        /// <param name="bases"></param>
        /// <returns></returns>
        public static float[] BaseFraction(IZmwBases bases)
        {
            return nChan.Fill(c => (float)-1.0);
        }
 
    }
}
