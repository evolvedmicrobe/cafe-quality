using System;
using Bio;
using Bio.IO.BAM;
using Bio.IO.SAM;
using System.Collections.Generic;
using System.Linq;
using Bio.IO.PacBio;
using System.Diagnostics;
using System.IO;
using VariantCaller;

namespace BAMErrorReporter
{
    class MainClass
    {
        static CCSWriter ccsw;
        static VariantWriter vwriter;
        static List<Variant> vempty = new List<Variant>();
        public static void Main (string[] args)
        {
            string[] bamNames = new string [] {
                "/Users/nigel/git/cafe-quality/NotTracked/C++-Zscore/master_26bbfc4_snr50000.bam",
                "/Users/nigel/git/cafe-quality/NotTracked/C++-Zscore/master_26bbfc4_snr2.bam",
                "/Users/nigel/git/cafe-quality/NotTracked/C++-Zscore/master_26bbfc4.bam"
                
            };
            foreach(var bamF in bamNames) {
                var fi = new FileInfo(bamF);
                var prefix = fi.Name.Split ('.') [0];
                var fname_prefix = Path.Combine(fi.DirectoryName, prefix);

                vwriter = new VariantWriter (fname_prefix + "_all_variants.csv");
                ccsw = new CCSWriter (fname_prefix + "_read_report.csv");
                var seqs = PacBioBamReader.ParseReads(bamF);
                foreach (var s in seqs) {
                    var r = ReferenceGenomes.AssignReadToReference (s);
                    if (r != null) {
                        var alns = r.AlignSequence (s.Sequence).ToList ();
                        if (alns.Count != 0) {
                            var bestS = alns.Max (z => z.Score);
                            var best = alns.Where (x => x.Score == bestS).First();
                            var variants = VariantCaller.VariantCaller.CallVariants (best, r.RefSeq, s.Sequence);
                            ccsw.Write (s, variants, r);
                            variants.ForEach(z => vwriter.Write(s, z));
                        }
                    }
                }
                ccsw.sw.Close ();
                vwriter.sw.Close ();
            }
        }
    }
}
