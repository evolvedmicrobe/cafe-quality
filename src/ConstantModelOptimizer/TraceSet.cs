//// TRANSLATED FROM EDNA
//using System;
//using System.Collections.Concurrent;
//using System.Collections.Generic;
//using System.IO;
//using System.Linq;
//using System.Threading;
//using System.Threading.Tasks;
//using PacBio.IO;
//using PacBio.Utils;
//using PacBio.Align;
//
//namespace ConstantModelOptimizer
//{
//    public class TraceSet
//    {
//
//        OldCmpH5Reader cmp;
//        BasCollection bas;
//        int maxTraces;
//
//        /// <summary>
//        /// A Basemap 
//        /// </summary>
//        Dictionary<char, int> bm;
//
//        Dictionary<string, float[]> blData;
//
//        Dictionary<string, float[]> sigmaData;
//
//        List<IAlnSummary> Alns;
//
//        public TraceSet(string cmpFile, string fofn, int maximumTraces)
//        {
//            cmp = OldCmpH5Reader.CreateReader (cmpFile);
//            bas = BasCollection.FromFofn (fofn);
//            maxTraces = maxTraces;
//            bm = new Dictionary<char, int> ();
//            for(int i=0; i< bas.BaseMap.Length; i++)
//            {
//                bm [i] = bas.BaseMap [i];
//            }
//            blData = readBaselineData (bm, fofn, out sigmaData);
//            var rng = new RandomCMWC(42);
//            Alns = cmp.Alns.TargetSample ((fun => true), maxTraces, cmp.Alns.Count, rng).ToList ();
//
//        }
//
//        public string ReportsFolder {
//            get { return cmp.ReportsFolder; }
//        }
//
//        public IAlignment GetAlignment(IAlnSummary aln)
//        {
//            return cmp.ReadAlignment (aln);
//        }
//
//        public float[] GetBaselineBias(string id)
//        {
//            blData [id];
//        }
//
//        public float[] GetBaselineSigma(string id)
//        {
//            return sigmaData [id];
//        }
//
//        public int Count {get{return Alns.Length;}}
//
//        public IEnumerable<Trace> Traces() 
//        {
//            foreach (var a in Alns) {
//                var zmwBases = bas.GetRead (a);
//                var tr = new Trace (this, zmwBases, a);
//                yield return tr;
//            }
//
//        }
//
//        public static void readBaseLineData(Dictionary<char, int> bm, string fofnFile, out Dictionary<string, float[]> sgData) 
//        {
//            var resultsDirs = System.IO.File.ReadAllLines (fofnFile).Select (x => Path.GetDirectoryName (x)).Distinct ();
//            var stsFiles = resultsDirs.SelectMany (d => System.IO.Directory.EnumerateFiles (d, "*.sts.csv"));
//
//            var blData = new Dictionary<string, float[]> ();
//            var sgData = new Dictionary<string, float[]>();
//
//            let readStsCsv f =
//                let csvRows = miniCsvReader f
//
//                for r in csvRows do
//                    let zmwId = (r("Movie")).Trim() + "/" + (r("Zmw")).Trim()
//
//                    let bias = Array.zeroCreate 4
//                    bias.[bm.['T']] <- float (r("BaselineLevel_T"))
//                    bias.[bm.['G']] <- float (r("BaselineLevel_G"))
//                    bias.[bm.['A']] <- float (r("BaselineLevel_A"))
//                    bias.[bm.['C']] <- float (r("BaselineLevel_C"))
//                    blData.[zmwId] <- bias
//
//                    let sigm = Array.zeroCreate 4
//                    sigm.[bm.['T']] <- float (r("BaselineStd_T"))
//                    sigm.[bm.['G']] <- float (r("BaselineStd_G"))
//                    sigm.[bm.['A']] <- float (r("BaselineStd_A"))
//                    sigm.[bm.['C']] <- float (r("BaselineStd_C"))
//                    sgData.[zmwId] <- sigm
//                    done
//
//
//                    stsFiles |> Seq.iter readStsCsv
//
//                    (blData, sgData)
//
//
//        }
//    }
//}
//
