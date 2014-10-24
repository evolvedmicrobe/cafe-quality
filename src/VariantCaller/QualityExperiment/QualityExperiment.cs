using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Threading.Tasks;

using Bio;
using Bio.IO.FastA;
using PacBio.IO;

namespace VariantCaller.QualityExperiment
{
    /// <summary>
    /// This class represents an experiment designed to evaluate the accuracy of CCS reads.  
    /// It consists of several pulse file and references.
    /// </summary>
    public class QualityExperiment
    {
        /// <summary>
        /// All references used in experiment
        /// </summary>
        public List<Reference> References;
        /// <summary>
        /// Original base files used.  
        /// </summary>
        public List<FileInfo> BaseFiles;

        /// <summary>
        /// A list of CCS Reads.
        /// </summary>
        public List<CCSRead> CCSReads;

        /// <summary>
        /// Make a new quality experiment
        /// </summary>
        /// <param name="baseFiles">The original bax files</param>
        /// <param name="ccsReads">The CCS consensus files</param>
        /// <param name="referenceFile">The refrence files</param>
        public QualityExperiment(IEnumerable<string> baseFiles, IEnumerable<string> ccsReads, IEnumerable<string> subReads,  string referenceFile)
        {
            // Verify Files exist
            var names = new[] {"baseFiles", "ccsReads", "subReads"};
            var files = new[] { baseFiles, ccsReads, subReads };
            for (int i=0; i<files.Length; i++)
            {
                var fss = files[i];
                if (fss != null)
                {
                    if (fss.Any(p => !(new FileInfo(p)).Exists))
                    {
                        throw new FileNotFoundException(names[i]);
                    }
                }
            }
            
            // Load Base files if present
            if (baseFiles != null)
            {
                BaseFiles = baseFiles.Select(x => new FileInfo(x)).ToList();
            }

            References = (new FastAParser(referenceFile)).Parse().Select(x => new Reference((Sequence)x)).ToList();

            // Start loading CCS reads
            var ccs_reads_maker = Task<List<CCSRead>>.Factory.StartNew( () =>
                {
                   return ccsReads.SelectMany(p => (new FastAZippedParser(p))
                                                .Parse()
                                                .Select(x => new CCSRead((Sequence)x)))
                            .ToList();
                });
            // Start loading the subreads
            var subReadDict = loadSubReads(subReads);
            
            // Merge
            CCSReads = ccs_reads_maker.Result;
			int missing = 0;
			List<CCSRead> missingReads = new List<CCSRead>();
            foreach (var read in CCSReads)
            {
				List<CCSSubRead> subReadList;
				bool hasValue = subReadDict.TryGetValue(read.ZMW, out subReadList);
				if (!hasValue) {
					missing++;
					missingReads.Add (read);
				} else {
					read.SubReads = subReadList;
				}
            }
			

            assignCCSReadsToReference ();

            var numMissing = CCSReads.Count ( x => x.AssignedReference == null);
            Console.WriteLine ("Total Reads: " + CCSReads.Count.ToString ());
            Console.WriteLine ("Unmatched between CCS and subreads (missing): " + missing.ToString () + " reads");
            var percMissing = numMissing / (double) CCSReads.Count;
            Console.WriteLine ("Not Assigned to References: " + numMissing.ToString () + " reads (" + percMissing.ToString("f4") +"%)");
            var counts = new int[References.Count];
            foreach (var read in CCSReads) {
                for(int j =0; j < counts.Length; j++) {
                    if (read.AssignedReference == References[j]) {
                        counts [j]++;
                    }
                }
            }
            for (int i = 0; i < References.Count; i++) {
                var c = counts [i];
                percMissing = c / (double) CCSReads.Count;
                Console.WriteLine ("Assigned to " + References[i].RefSeq.ID + ": " + c.ToString () + " reads (" + percMissing.ToString("f4") +"%)");

            }

        }
        /// <summary>
        /// Naive temporary implementation, where we assign on size.
        /// TODO: Improvements needed!
        /// </summary>
        private void assignCCSReadsToReference()
        {
            foreach (var v in CCSReads) {
                foreach (var r in References) {
                    if (Math.Abs (r.RefSeq.Count - v.Seq.Count) < 25) {
                        v.AssignedReference = r;
                        break;
                    }
                }
            }
        }

        /// <summary>
        /// Load the subreads in parallel.
        /// </summary>
        /// <param name="subReads"></param>
        /// <returns></returns>
        private SubReadCollection loadSubReads(IEnumerable<string> subReads)
        {
            var subReadDictionaries = subReads.AsParallel().Select(z => new SubReadCollection(z)).ToList();
            var mainDictionary = subReadDictionaries[0];
            for (int i = 1; i < subReadDictionaries.Count; i++)
            {
                var curDict = subReadDictionaries[i];
                foreach (var kv in curDict) {
                    mainDictionary[kv.Key] = kv.Value;
                }
                subReadDictionaries[i] = null;
            }
            return mainDictionary;
        }

        public IEnumerable<IZmwBases> GetZMWReads()
        {
            foreach (var bas in BaseFiles)
            {

                var br = BaseReader.CreateSource(bas.FullName);
                // null range apparently selects all
                foreach (var zmw in br.ByHoleNumberRange(null))
                {
                    yield return zmw;
                }
            }
        }

    }    
 }
