using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Threading.Tasks;

using Bio;
using Bio.IO.FastA;
using PacBio.Data;
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
            foreach (var read in CCSReads)
            {
                read.SubReads = subReadDict[read.ZMW];
            }
            subReadDict = null;

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
