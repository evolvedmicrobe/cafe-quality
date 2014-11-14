using System;
using VariantCaller;
using Bio;
using System.Linq;
using System.Collections.Generic;
using System.Diagnostics;

namespace HPcorrector
{
    [DebuggerDisplay("Length = {Length}, BP = {BP}")]
    public struct Homopolymer {
        public int Start;
        public int Length;
        public char BP;
        public Homopolymer(int start, int length, byte bp)
        {
            this.Start = start;
            this.Length = length;
            this.BP = (char)bp;
        }
    }

    public class CCSReadCorrector
    {
        private const int MINIMUM_HP_LENGTH = 4;
        public const int MINIMUM_READS_NEEDED = 4;
        public const double MINIMUM_ERROR_PERCENTAGE = 0.1;
        public static int Fix_Count = 0;
        public const double MIN_ACCEPTABLE_RATIO = .5;
        public static CCSRead CorrectRead(CCSRead read)
        {
            var hps = CallHomopolymers (read);
            if (hps.Count > 0) {
                hps.Reverse ();
                var orgSeq = read.Seq.GetInternalArray().ToList();
                foreach (var hp in hps) {
                   
                    var needsFix = DecideIfHPNeedsFixByRatio (read, hp);
                    if (needsFix) {
                        System.Threading.Interlocked.Increment (ref Fix_Count);
                        if (Fix_Count % 100 == 0) {
                            Console.WriteLine ("Fixed: " + Fix_Count.ToString ());
                        }
                        orgSeq.Insert (hp.Start, (byte)hp.BP);
                    }
                }
                var orgID = read.Seq.ID;
                read.Seq = new Sequence (DnaAlphabet.Instance, orgSeq.ToArray(), false);
                read.Seq.ID = orgID;
            }
            return read;
        }

        private static List<Homopolymer> CallHomopolymers(CCSRead read)
        {
            var seq = read.Seq.GetInternalArray ();
            var lastBase = seq[0];
            var homos = new List<Homopolymer> ();
            for (int i = 1; i < seq.Length; i++) {
                if (seq [i] == lastBase) {
                    int start = i - 1;
                    int len = 2;
                    while (i < (seq.Length - 1) && seq [i + 1] == lastBase) {
                        len += 1;
                        i++;                        
                    }
                    if (len >= MINIMUM_HP_LENGTH) {
                        var hp = new Homopolymer (start, len, lastBase);
                        homos.Add (hp);
                    }
                } else {
                    lastBase = seq [i];
                }
            }
            return homos;
        }
       // private static System.IO.StreamWriter sw = new System.IO.StreamWriter("/Users/nigel/git/cafe-quality/data/ratios.csv");
        private static bool DecideIfHPNeedsFixByRatio(CCSRead read, Homopolymer hp)
        {
            //var cutPoint = hp.BP == (byte)'A' || hp.BP == (byte)'T' ? .333 : .5;
            // First to align all subreads and make a decision
            var alner = new Reference (read.Seq);
            var aln2 = new Reference (read.AssignedReference.RefSeq.GetReverseComplementedSequence () as Sequence);
            bool result = false;
            //TODO: Obviously repeating this for all subreads is insane
            int delCount = 0;
            int insCount = 0;
            int alignments = 0;
            double rat = 0;
            foreach (var subRead in read.SubReads) {
                if ( Math.Abs(subRead.Seq.Length - read.Seq.Count) > 25) {
                    continue;
                }
                var aln =  alner.AlignSequence(new Bio.Sequence(DnaAlphabet.Instance, subRead.Seq, false));
                if (aln.Count > 0) {
                    Console.WriteLine ("To Seq");
                    Console.WriteLine(aln[0].ToString());
                    Console.WriteLine ("To Ref");
                    Console.WriteLine(aln2.AlignSequence(new Bio.Sequence(DnaAlphabet.Instance, subRead.Seq, false))[0].ToString());

                    alignments++;
                    var variants = VariantCaller.VariantCaller.CallVariants (aln [0], read.AssignedReference.RefSeq).Where(x => x.StartPosition == hp.Start -1).ToList();
                    if (variants.Count > 0 && variants [0].Type == VariantType.INDEL) {
                        var type = variants [0] as IndelVariant;
                        Console.WriteLine (hp.BP.ToString () + " - " + hp.Length.ToString ());
                        if (type.InsertionOrDeletion == IndelType.Insertion) {
                            Console.WriteLine ("Insertion");
                            Console.WriteLine (aln [0].ToString ());
                            insCount++;
                        } else {
                            Console.WriteLine ("Deletion");
                            Console.WriteLine (aln [0].ToString ());
                            delCount++;
                        }
                    } 

                }
            }
            var totalErrors = ((double)(insCount + delCount));
            rat = (double)delCount / totalErrors;
            var percErrors = totalErrors /  alignments;
            if (delCount + insCount >= MINIMUM_READS_NEEDED) {
                result = rat <= MIN_ACCEPTABLE_RATIO && percErrors > MINIMUM_ERROR_PERCENTAGE;
            }
          //  sw.WriteLine (read.AssignedReference.RefSeq.ID+","+read.OriginallyRevComped+","+read.ZMWnumber+","+hp.Start+","+hp.Length.ToString () + "," + hp.BP.ToString () + "," + rat.ToString ()+","+insCount.ToString()+","+delCount.ToString()+","+read.SubReads.Count+","+alignments);
           // sw.Flush ();


            return result;
        }
        private static bool DecideIfHPNeedsFixByVoting(CCSRead read, Homopolymer hp)
        {
            //var cutPoint = hp.BP == (byte)'A' || hp.BP == (byte)'T' ? .333 : .5;
            // First to align all subreads and make a decision
            var alner = new Reference (read.Seq);
            bool result = false;
            //TODO: Obviously repeating this for all subreads is insane
            int delCount = 0;
            int insCount = 0;
            int normCount = 0;
            int alignments = 0;
            double rat = 0;
            foreach (var subRead in read.SubReads) {
                if (Math.Abs (subRead.Seq.Length - read.Seq.Count) > 25) {
                    continue;
                }
                var aln = alner.AlignSequence (new Bio.Sequence (DnaAlphabet.Instance, subRead.Seq, false));
                if (aln.Count > 0) {
                    alignments++;
                    var variants = VariantCaller.VariantCaller.CallVariants (aln [0], read.Seq).Where (x => x.StartPosition == hp.Start - 1).ToList ();
                    if (variants.Count == 0 && aln [0].FindQueryPositionCorrespondingtoReferencePosition (hp.Start - 1).HasValue) {
                        normCount++;
                    } else if (variants.Count > 0 && variants [0].Type == VariantType.INDEL) {
                        var type = variants [0] as IndelVariant;
                        if (type.InsertionOrDeletion == IndelType.Insertion) {
                            insCount++;
                        } else {
                            delCount++;
                        }
                    }
                }
            }
            if (alignments == 0) {
                return false;
            }
            var totalCov = ((double)(insCount + delCount + normCount));
            rat = (double)insCount / totalCov;
            return rat > 0.5;

        }

        public static void CountAlns(CCSRead read)
        {
            var alner = new Reference (read.AssignedReference.RefSeq);
            int sub = 0;
            Console.WriteLine (read.SubReads.Count);
            foreach (var subRead in read.SubReads) {
                sub++;
                if ( Math.Abs(subRead.Seq.Length - read.Seq.Count) > 25) {

                    Console.WriteLine (sub.ToString () + ",LengthDiff");
                    continue;
                }
                var aln =  alner.AlignSequence(new Bio.Sequence(DnaAlphabet.Instance, subRead.Seq, false));
                if (aln.Count > 0) {
                    //Console.WriteLine (aln [0].ToString ());
                    Console.WriteLine (sub.ToString () + ",Exists");

                } else {
                    Console.WriteLine (sub.ToString () + ",Null");
                }
            }
           
        }
    }
}

