using System;
using System.Linq;
using NUnit.Framework;
using PacBio.Utils;

namespace PacBio.Consensus.Test
{
    [TestFixture]
    public class POATest
    {
        [Test]
        public void ConsensusCorePoaTest()
        {
            var s1 = "TTTACAGGATAGTCCAGT";
            var s2 = "ACAGGATACCCCGTCCAGT";
            var s3 = "ACAGGATAGTCCAGT";
            var s4 = "TTTACAGGATAGTCCAGTCCCC";
            var s5 = "TTTACAGGATTAGTCCAGT";
            var s6 = "TTTACAGGATTAGGTCCCAGT";
            var s7 = "TTTACAGGATAGTCCAGT";

            var reads = new []{s1, s2, s3, s4, s5, s6, s7};

            float score;
            var consensus = ConsensusCorePoa.FindConsensus(reads, out score);


            Assert.AreEqual("TTTACAGGATAGTCCAGT", consensus);
        }

        [Test]
        public void LocalStaggered()
        {
            var s1 = "TTTACAGGATAGTGCCGCCAATCTTCCAGT";
            var s2 =    "GATACCCCGTGCCGCCAATCTTCCAGTATATACAGCACGAGTAGC";
            var s3 =         "ATAGTGCCGCCAATCTTCCAGTATATACAGCACGGAGTAGCATCACGTACGTACGTCTACACGTAATT";
            var s4 =                                                             "ACGTCTACACGTAATTTTGGAGAGCCCTCTCTCACG";
            var s5 =                                                                   "ACACGTAATTTTGGAGAGCCCTCTCTTCACG";
            var s6 =      "AGGATAGTGCCGCCAATCTTCCAGTAATATACAGCACGGAGTAGCATCACGTACG";
            var s7 =         "ATAGTGCCGCCAATCTTCCAGTATATACAGCACGGAGTAGCATCACGTACGTACGTCTACACGT";
            var res=         "ATAGTGCCGCCAATCTTCCAGTATATACAGCACGGAGTAGCATCACGTACGTACGTCTACACGTAATT".Replace(" ", "");

            var reads = new string[] {s1, s2, s3, s4, s5, s6, s7 };

            var poa = new PoaLocal();
            poa.AddReads(reads);

            float consensusScore;
            string consensus;
            poa.FindConsensusAndAlignments(1, out consensusScore, out consensus);
            Console.WriteLine(consensus);

            Assert.AreEqual(res, consensus);
        }


        [Test]
        public void LocalLongInsert()
        {
            var s1 = "TTTACAGGATAGTGCCGCCAATCTTCCAGTGATACCCCGTGCCGCCAATCTTCCAGTATATACAGCACGAGGTAGC";
            var s2 = "TTTACAGGATAGTGCCGGCCAATCTTCCAGTGATACCCCGTGCCGCCAATCTTCCAGTATATACAGCACGAGTAGC";
            var s3 = "TTGTACAGGATAGTGCCGCCAATCTTCCAGTGATGGGGGGGGGGGGGGGGGGGGGGGGGGGACCCCGTGCCGCCAATCTTCCAGTATATACAGCACGAGTAGC";

            var reads = new string[] { s1, s2, s3 };

            var poa = new PoaLocal();
            poa.AddReads(reads);

            float consensusScore;
            string consensus;
            poa.FindConsensusAndAlignments(1, out consensusScore, out consensus);
            Console.WriteLine(consensus);

            Assert.AreEqual("TTTACAGGATAGTGCCGCCAATCTTCCAGTGATACCCCGTGCCGCCAATCTTCCAGTATATACAGCACGAGTAGC", consensus);
        }
    }   
}
