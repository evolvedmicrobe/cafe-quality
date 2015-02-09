using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;

namespace PacBio.Consensus
{
    /// <summary>
    /// Crappy hard coding for this testing
    /// </summary>
    public static class WhiteList
    {
        static string fname = @"/Users/nigel/16S/BEI/p6c4_BEI.fixed.zmws";
        static HashSet<string> okayZMWs = new HashSet<string>();

        static WhiteList ()
        {
            var lines = File.ReadLines (fname).Select (z => String.Join ("/", z.Split('/').Take (2).ToArray()));
            foreach (var l in lines) {
                okayZMWs.Add (l);
            }

        }
        public static bool ZMWisOkay(string movieAndHole)
        {
            return okayZMWs.Contains (movieAndHole);
        }
    }
}

