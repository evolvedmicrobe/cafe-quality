using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using Bio;

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
            var name = seq.ID;
			var Seq = seq.ToArray ();
            var sp = name.Split('/');
            ParentZMW = Convert.ToInt32(sp[1]);
            var spl = sp.Last();
            var spls = spl.Split(' ');
            RQ = Convert.ToSingle(spls.Last().Split('=').Last());
            var pos = spls.First().Split('_').Select(z => Convert.ToUInt32(z)).ToArray();
            Start = pos[0];
            End = pos[1];
        }
    }
}
