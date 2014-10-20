using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using Bio;
using Bio.IO.FastA;

namespace VariantCaller.QualityExperiment
{
    /// <summary>
    /// A dictionary that maps ZMWs to a list of their subreads.
    /// </summary>
    public class SubReadCollection : IEnumerable<KeyValuePair<int,List<CCSSubRead>>>
    {
        Dictionary<int, List<CCSSubRead>> zmwsToSubReads = new Dictionary<int, List<CCSSubRead>>(10000);
        public SubReadCollection(string subReadFile)
        {

            var reads = (new FastAZippedParser(subReadFile).Parse())
                .Select(p => new CCSSubRead((Sequence)p));
            List<CCSSubRead> curList = null;
            var lastZMW = Int32.MinValue;
            foreach(var subRead in reads)
            {
                var zmw = subRead.ParentZMW;
                if (zmw == lastZMW)
                {
                    curList.Add(subRead);
                }
                else
                {
                    // Should never happen if reads are in order
                    var alreadyThere = zmwsToSubReads.TryGetValue(zmw, out curList);
                    if (!alreadyThere)
                    {
                        curList = new List<CCSSubRead>();
                        zmwsToSubReads[zmw] = curList;
                    }
                    curList.Add(subRead);
                }
                lastZMW = zmw;
            }
        }

        /// <summary>
        /// Get a list of subreads for a given ZMW.
        /// </summary>
        /// <param name="zmw"></param>
        /// <returns></returns>
        public List<CCSSubRead> this[int zmw] {
            get { return zmwsToSubReads[zmw]; }
            internal set { zmwsToSubReads[zmw] = value; }
        }

        public IEnumerator<KeyValuePair<int, List<CCSSubRead>>> GetEnumerator() {
            return zmwsToSubReads.GetEnumerator();
        }

		public bool TryGetValue(int zmw, out List<CCSSubRead> values)
		{
			return zmwsToSubReads.TryGetValue (zmw, out values);
		}
        System.Collections.IEnumerator System.Collections.IEnumerable.GetEnumerator()
        {
            return this.GetEnumerator();
        }
    }
}
