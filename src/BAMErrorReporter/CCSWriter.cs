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
    public class CCSWriter
    {
        public StreamWriter sw;
        public CCSWriter (string filename)
        {
            sw = new StreamWriter (filename);
            sw.WriteLine ("Movie,ZMW,Reference,NumSubReads,NumErrors,Length,NumIndelErrors,NumSNPErrors,SnrT,SnrG,SnrA,SnrC,Zg,Za,Success,MemFail,AlphaBeta,LowZ,RQ,ComputeMS,EndVariants");

        }

        public void Write(PacBioCCSRead read, List<Variant> variants, Reference assignedRef) {
            var reference = assignedRef.RefSeq.ID;
            var indel_cnt = variants.Count (x => x.Type == VariantType.INDEL);
            var snp_cnt = variants.Count (x => x.Type == VariantType.SNP);
            var end_var = variants.Count (x => x.AtEndOfAlignment);
            var tow = String.Join (",", 
                          read.Movie,
                          read.HoleNumber.ToString (),
                          reference,
                          read.NumPasses.ToString (),
                          variants.Count.ToString (),
                          read.Sequence.Count.ToString (),
                          indel_cnt.ToString (),
                          snp_cnt.ToString (),
                          read.SnrT.ToString (),
                          read.SnrG.ToString (),
                          read.SnrA.ToString (),
                          read.SnrC.ToString (),
                          read.GlobalZscore.ToString (),
                          read.AvgZscore.ToString (),
                          read.ReadsSuccessfullyAdded.ToString (),
                          read.ReadsMemFail.ToString(),
                          read.ReadsAlphaBetaMismatch.ToString (),
                          read.ReadsBadZscore.ToString (),
                          read.ReadQuality.ToString (),
                          read.ComputingMilliSeconds,
                          end_var.ToString ());
            sw.WriteLine (tow);
            sw.Flush ();
        }

        private string homopolymerLength(Variant v)
        {
            if (v.AtEndOfAlignment)
                return "-999";
            if (v is IndelVariant) {
                return (v as IndelVariant).HomopolymerLengthInReference.ToString ();
            } else if (v is SNPVariant) {
                return "1";
            }
            throw new Exception ("Unknown variant type");

        }

    }
}

