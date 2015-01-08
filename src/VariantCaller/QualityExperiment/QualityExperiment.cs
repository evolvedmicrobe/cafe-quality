using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Threading.Tasks;

using Bio;
using Bio.IO.FastA;
using PacBio.Data;

namespace VariantCaller
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
        public List<FileInfo> BaseFileNames;

        private List<BaxH5Reader> baseFileReaders;

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
                BaseFileNames = baseFiles.Select(x => new FileInfo(x)).ToList();
                baseFileReaders = BaseFileNames.Select (x => new PacBio.Data.BaxH5Reader(x.FullName)).ToList ();
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
				bool hasValue = subReadDict.TryGetValue(read.ZMWnumber, out subReadList);
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
            CCSRead.ParentExperiment = this;
        }
        /// <summary>
        /// Naive temporary implementation, where we assign on size.
        /// TODO: Improvements needed!
        /// </summary>
        private void assignCCSReadsToReference()
        {
            var lambda = References.Find(r => r.RefSeq.ID.StartsWith("lambda_NEB3011", StringComparison.Ordinal));
            foreach (var v in CCSReads) {
                if (v.Movie.StartsWith ("m141115", StringComparison.Ordinal)) {
                    if (v.Seq.Count > 60) {
                        v.AssignedReference = lambda;
                    }
                } else {
                    foreach (var r in References) {
                        if (Math.Abs (r.RefSeq.Count - v.Seq.Count) < 25) {
                            v.AssignedReference = r;
                            break;
                        }
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

        public IEnumerable<Zmw> GetZMWSequencingReads()
        {
            foreach (var bas in baseFileReaders)
            {
                foreach (var zmw in bas.SequencingZmws)
                {
                    yield return zmw;
                }
            }
        }

        /// <summary>
        /// Because subreads are read from a fasta file produced by consensus tools, which
        /// uses PacBio.IO, whilst we are now using PacBio.Data, they may not match.  Particularly
        /// since the adapter finding step is unique to Pat's code.
        /// 
        /// This procedure matches the subreads from the fasta file with those from the 
        /// PacBio.Data class, it makes sure the sequences match.
        /// </summary>
        /// <param name="read">Read.</param>
        public  void ValidateSubReads(CCSRead read)
        {
            if (read.ZMW == null) {
                read.ZMW = CCSRead.ParentExperiment.GetZMWforRead(read);
            }

        }

        public Zmw GetZMWforRead(CCSRead read)
        {

            foreach (var v in baseFileReaders) {
                if (v.HasHoleNumber (read.ZMWnumber)) {
                    return v [read.ZMWnumber];
                }
            }
            throw new BioinformaticsException ("The loaded experiment did not have data for that hole");
        }

    }    
 }
