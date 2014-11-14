using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
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
}
