using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using Bio;
using Bio.Util;

namespace VariantCaller
{
    public class CCSSubRead
    {
        public readonly uint Start, End;
        public readonly byte[] Seq;
        public readonly float RQ;
        public readonly int ParentZMW;

        public CCSSubRead(Sequence seq)
        {
			Seq = seq.GetInternalArray();
			// >m141008_060349_42194_c100704972550000001823137703241586_s1_p0/54494/0_116 RQ=0.866
            var name = seq.ID;
			var useful = FastStringUtils.ReturnStringAfterDelimiter (name, '/');
			string left, right;
			FastStringUtils.SplitStringIntoTwo (useful, '/', out left, out right);
			ParentZMW = Convert.ToInt32(left);
            
			// Now to split the last bit (e.g. 0_116 RQ=0.866)
			int pos_underscore = 0;
			int pos_space = 0;
			int pos_equals = 0;
			for (int i = 0; i < right.Length; i++) {
				var curC = right [i];
				if (curC == '_') {
					pos_underscore = i;
				}
				else if (curC == ' ') {
					pos_space = i;
				}
				else if (curC == '=') {
					pos_equals = i;
					break;
				}
			}

			var start_s = right.Substring (0, pos_underscore);
			var end_s = right.Substring (pos_underscore+1, pos_space - pos_underscore - 1);
			var rq_s = right.Substring (pos_equals + 1, right.Length - pos_equals - 1);

            RQ = Convert.ToSingle(rq_s);
			Start = Convert.ToUInt32(start_s);
			End = Convert.ToUInt32(end_s);
        }
    }
}
