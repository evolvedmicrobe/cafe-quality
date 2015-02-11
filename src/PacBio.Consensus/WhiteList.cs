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

        /// <summary>
        /// These are species in the training data that have several similar sequences present, because of this we will only train using data with identical sequences, as shown below.
        /// 
        /// See the python script /Users/nigel/16S/ExamineFasta.py for more details on this issue.
        /// </summary>
        static HashSet<string> okTemplates = new HashSet<string>() { "ALL4MER.V2.01"};

        public static bool TemplateIsOkay(string name)
        {

            return okTemplates.Contains (name);
        }
    }
}