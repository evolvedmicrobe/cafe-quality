using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using NUnit.Framework;
using PacBio.Align;
using PacBio.IO.Fasta;
using PacBio.Utils;

namespace PacBio.Consensus
{
    /// <summary>
    /// Code in this file is for generating sets of barcodes that are well spaced from each other, and optionally for the RC's of each other.
    /// Used to generate the new set of 450 barcodes.
    /// </summary>
    internal class GenerateBarcodes
    {
        /// <summary>
        /// You can often find many more barcodes by cycling through the random seed, especially for short barcodes. 
        /// </summary>
        public GenerateBarcodes(int randomSeed, int barcodeLength)
        {
            r = new RandomCMWC(randomSeed);
            this.length = barcodeLength;
        }

        private int length = -1;
        private RandomCMWC r = new RandomCMWC(0);
        private char[] bases = new char[] {'A', 'C', 'G', 'T'};

        /// <summary>
        /// Generate random barcodes with max HP Length of 2
        /// </summary>
        public string MakeCandidate()
        {
            var c = new char[length];

            for (int i = 0; i < length; i++)
            {
                var newBase = bases[r.Next(4)];

                while (i > 0 && newBase == c[i - 1])
                //while (i > 1 && c[i - 1] == c[i - 2] && newBase == c[i - 1])
                    newBase = bases[r.Next(4)];

                c[i] = newBase;
            }

            return new string(c);
        }

        // Edit distance scores
        // Indels are counted as 2, mismatches as 3 -- favors mismatches a bit because we are less likely to make that 
        // error
        public readonly short[] scores = new short[] {-2, -2, -3, 0};

        public int GetEditDistance(string s1, string s2)
        {
            return -GlobalAlign.GetGlobalAlignScore(s1, s2, scores);
        }

        public bool ScreenCandidate(string candidate, List<string> barcodes, int minEditDist, bool checkReverseComplement)
        {
            var rcCandidate = DNA.ReverseComplement(candidate);

            var hasExistingCloseMatch = barcodes.AsParallel().Any(
                existing =>
                    {
                        var ed = GetEditDistance(candidate, existing);

                        if (ed < minEditDist)
                        {
                            return true;
                        }

                        if (checkReverseComplement)
                        {
                            ed = GetEditDistance(rcCandidate, existing);

                            if (ed < minEditDist)
                            {
                                return true;
                            }
                        }

                        return false;
                    });

            return !hasExistingCloseMatch;
        }

        public List<string> DoDesign(int nbarcodes, int startEditDistance, int finalEditDistance, bool checkReverseComplement)
        {
            var list = new List<string>();

            list.Add(MakeCandidate());
            var minEditDistance = startEditDistance;
            var maxTries = 300000;

            while (list.Count < nbarcodes)
            {
                var c = MakeCandidate();
                var nTried = 1;

                while (!ScreenCandidate(c, list, minEditDistance, checkReverseComplement))
                {
                    c = MakeCandidate();
                    nTried++;

                    if (nTried > maxTries)
                    {
                        break;
                    }
                }

                // Only add it if it's GTG
                if (ScreenCandidate(c, list, minEditDistance, checkReverseComplement))
                {
                    list.Add(c);
                    Console.WriteLine("Got Barcode {0}, MinEd: {1} -- Tries: {2}", list.Count, minEditDistance, nTried);
                }

                if (nTried > maxTries && minEditDistance == finalEditDistance)
                {
                    Console.WriteLine("Couldn't find new barcodes at minEditDistance, exiting");
                    break;
                }

                if (nTried > 0.75 * maxTries && minEditDistance > finalEditDistance)
                {
                    Console.WriteLine("Having trouble finding sequences. Dropping minEditDistance");
                    minEditDistance--;
                }
            }

            return list;
        }
    }

    [TestFixture, Explicit]
    public class BarcodeDesignTest
    {
        /// <summary>
        /// Generates a set of ~450 16bp barcodes, disallowing barcodes that are RC of each other
        /// </summary>
        [Test]
        public void GenerateBigBarcodeSetNoRC()
        {
            
            var bcd = new GenerateBarcodes(12, 16);
            var barcodes = bcd.DoDesign(450, 16, 13, true);

            var lines = new List<string>();
            int count = 1;

            foreach (var bc in barcodes)
            {
                lines.Add(String.Format(">lbc{0}", count));
                lines.Add(bc);
                count++;
            }

            System.IO.File.WriteAllLines(@"barcodes-big.fasta", lines.ToArray());
        }


        /// <summary>
        /// Generates a set of ~50 7bp barcodes, allowing barcodes that are RC of each other.
        /// For use as 'part of the adapter' barcodes.
        /// </summary>
        [Test]
        public void GenerateSmallBarcodeSetAllowRC()
        {
            var bcd = new GenerateBarcodes(1, 8);
            var barcodes = bcd.DoDesign(50, 11, 8, false);

            var lines = new List<string>();
            int count = 1;

            foreach (var bc in barcodes)
            {
                lines.Add(String.Format(">shortbc{0}", count));
                lines.Add(bc);
                count++;
            }

            System.IO.File.WriteAllLines(@"adapter-barcodes.fasta", lines.ToArray());
        }


        [Test]
        public void CheckEditDistance()
        {
            var bcd = new GenerateBarcodes(0, 16);

            var s1 = "ACGTACGTACGTACGT";
            var s2 = "ACGTCCGTATGTACGT";

            var ed1 = bcd.GetEditDistance(s1, s2);
            Assert.AreEqual(2, ed1);

            s1 = "ACGTACGTACGTACGT";
            s2 = "ACGTCGTCCGTACGGT";

            ed1 = bcd.GetEditDistance(s1, s2);
            Assert.AreEqual(3, ed1);


            s1 = "ACGTACGTACGTACGT";
            s2 = "ACGTTTTTACGTACGT";

            ed1 = bcd.GetEditDistance(s1, s2);
            Assert.AreEqual(3, ed1);


            s1 = "ACGTACGTACGTACGT";
            s2 = "ACTACGCTACCGTACT";

            ed1 = bcd.GetEditDistance(s1, s2);
            Assert.AreEqual(4, ed1);

            s1 = "ACGTACGTACGTACGT";
            s2 = "ACGTACGTACGT";

            ed1 = bcd.GetEditDistance(s1, s2);
            Assert.AreEqual(4, ed1);


            s1 = "ACACACATATGTGTACAGATACAGT";
            s2 = "ACACCCATATGTGTACAGTTACCCCAGT";

            ed1 = bcd.GetEditDistance(s1, s2);
            Assert.AreEqual(5, ed1);

        }

        /// <summary>
        /// Double check that the barcodes actually meet 
        /// </summary>
        [Test]
        public void ValidateBarcodes()
        {
            // original barcodes
            //var bcFile = @"\\usmp-acid\pmarks\Devel\Barcodes\barcode_complete.fasta";

            // new barcodes
            //var bcFile = @"c:\users\pmarks\barcodes-nohp-ed6.fasta";
            //var bcFile = @"c:\users\pmarks\barcodes-t1.fasta";
            var bcFile = @"C:\Users\pmarks\barcodes-ed65-450.fasta";


            var r = new SimpleFASTAReader(bcFile);
            var barcodes = r.Select(s => s.GetSequence()).ToArray();
            var bcd = new GenerateBarcodes(1, 16);

            var bestScore = int.MinValue;

            // Change this if you made barcodes w/o RC screening.
            var checkReverseComplement = true;


            for (int i = 0; i < barcodes.Length; i++)
            {
                //var score = GlobalAlign.GetGlobalAlignScore(barcodes[i], barcodes[i].ReverseComplement(), bcd.scores);


                for (int j = 0; j < barcodes.Length; j++)
                {
                    if (i == j)
                        continue;

                    var score = GlobalAlign.GetGlobalAlignScore(barcodes[i], barcodes[j], bcd.scores);

                    if (score > bestScore)
                        bestScore = score;

                    if (checkReverseComplement)
                    {
                        score = GlobalAlign.GetGlobalAlignScore(DNA.ReverseComplement(barcodes[i]), barcodes[j],
                                                                bcd.scores);

                        if (score > bestScore)
                            bestScore = score;
                    }
                }
            }

            Console.WriteLine("Checked all pairs of {0} barcodes.  Max Global Alignment Score: {1}", barcodes.Length, bestScore);
        }


        [Test]
        public void FindDups()
        {
            var bcFile1 = @"C:\Users\pmarks\barcodes-ed65-450.fasta";
            var bc2 = @"\\usmp-acid\pmarks\Devel\Barcodes\barcode_complete.fasta";

            var s1 = (new SimpleFASTAReader(bcFile1)).Select(r => r.GetSequence()).ToArray();
            var s2 = (new SimpleFASTAReader(bc2)).Select(r => r.GetSequence()).ToArray();

            var set1 = new HashSet<string>(s1);
            var set2 = new HashSet<string>(s2);

            
            set1.IntersectWith(set2);

            if (set1.Count > 0)
            {
                Console.WriteLine("Has overalps");
            }
        }
    }
}
