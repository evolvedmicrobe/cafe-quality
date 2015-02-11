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
        public static StreamWriter sw = new StreamWriter(@"/Users/nigel/16S/BEI/species.csv");

        static HashSet<string> okayZMWs = new HashSet<string>();

        /// <summary>
        /// These are species in the training data that have several similar sequences present, because of this we will only train using data with identical sequences, as shown below.
        /// 
        /// See the python script /Users/nigel/16S/ExamineFasta.py for more details on this issue.
        /// </summary>
        static HashSet<string> okSpecies = new HashSet<string>() {"N.meningitidis", "H.pylori", "A.baumannii", "P.gingivalis", "L.gasseri", "S.pneumoniae", "A.odontolyticus", "S.agalactiae"};

        static WhiteList ()
        {
            sw.WriteLine ("Species");
            var lines = File.ReadLines (fname).Select (z => String.Join ("/", z.Split('/').Take (2).ToArray()));
            foreach (var l in lines) {
                okayZMWs.Add (l);
            }

        }
        public static bool ZMWisOkay(string movieAndHole)
        {
            return okayZMWs.Contains (movieAndHole);
        }
        public static bool SpeciesIsOkay(string speciesName)
        {
            var newName = String.Join (".", speciesName.Split ('.').Take (2).ToArray ());
            lock (sw) {
                sw.WriteLine (newName);
            }
            return okSpecies.Contains (newName);
        }
    }
}

