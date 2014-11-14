using System;
using System.IO;
using System.Linq;
using VariantCaller;
using System.Threading.Tasks;
using System.Threading;
namespace HPcorrector
{
    class MainClass
    {
        public static void Main (string[] args)
        {
            var exp = GetExperiment ();
            var newOutFile = @"/Users/nigel/git/cafe-quality/data/corrected_ccs_ratioRef0.5.fa";
            var fastaOut = new Bio.IO.FastA.FastAFormatter (newOutFile);
            int count = 0;

            Parallel.ForEach (exp.CCSReads, ccs => {
            //foreach (var ccs in exp.CCSReads) {
                //CCSReadCorrector.CountAlns(ccs);
               // return;
                if (ccs.AssignedReference != null && ccs.AssignedReference.RefSeq.ID != "SmrtBellSequence") {
                    var alns = ccs.AssignedReference.AlignSequence (ccs.Seq);
                    if (alns.Count == 0) {
                        ccs.OriginallyRevComped = "NoAln";
                    } else {
                        var top = alns [0];
                        if (top.SecondSequence.Metadata.ContainsKey ("+isReversed")) {
                            ccs.OriginallyRevComped = "Reversed";
                        } else {
                            ccs.OriginallyRevComped = "Forward";
                        }
                    }
                
                    var newCCS = CCSReadCorrector.CorrectRead (ccs);
                    if (String.IsNullOrWhiteSpace (ccs.Seq.ID)) {
                            throw new Exception ();
                        }
                        lock (fastaOut) {
                            Interlocked.Increment (ref count);
                            if (count % 100 == 0) {
                                Console.WriteLine (count.ToString ());
                            }
                            fastaOut.Write (newCCS.Seq);
                        }
                    }
                });
                fastaOut.Close ();
            }
        public static QualityExperiment GetExperiment()
        {
            var direc = @"/Users/nigel/git/cafe-quality/data/";
            var ccsFiles = (new DirectoryInfo (direc)).GetFiles ().Where (h => h.Name.EndsWith (".ccs.fasta.gz"))
                .Select (u => u.FullName).ToList ();


            var subReads = (new DirectoryInfo (direc)).GetFiles ().Where (h => h.Name.EndsWith (".subreads.fasta.gz")) 
                .Select (u => u.FullName).ToList ();

            var reference = Path.Combine (direc, "References.fna");

            var qc_exp = new QualityExperiment (null, ccsFiles, subReads, reference);
            return qc_exp;    
        }
    }
}
