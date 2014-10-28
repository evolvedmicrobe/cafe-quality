using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Bio;
using VariantCaller;


namespace TestAlgorithms
{
    [TestClass]
    public class VariantCallTests
    {
        public readonly static string rCRS =
            @"AGTAGGCTCCCTTCCCCTACTCATCGCACTAATTTACACTCACAACACCCTAGGCTCACT
AAACATTCTACTACTCACTCTCACTGCCCAAGAACTATCAAACTCCTGAGCCAACAACTT
AATATGACTAGCTTACACAATAGCTTTTATAGTAAAGATACCTCTTTACGGACTCCACTT
ATGACTCCCTAAAGCCCATGTCGAAGCCCCCATCGCTGGGTCAATAGTACTTGCCGCAGT
ACTCTTAAAACTAGGCGGCTATGGTATAATACGCCTCACACTCATTCTCAACCCCCTGAC
AAAACACATAGCCTACCCCTTCCTTGTACTATCCCTATGAGGCATAATTATAACAAGCTC
CATCTGCCTACGACAAACAGACCTAAAATCGCTCATTGCATACTCTTCAATCAGCCACAT
AGCCCTCGTAGTAACAGCCATTCTCATCCAAACCCCCTGAAGCTTCACCGGCGCAGTCAT
TCTCATAATCGCCCACGGGCTTACATCCTCATTACTATTCTGCCTAGCAAACTCAAACTA
CGAACGCACTCACAGTCGCATCATAATCCTCTCTCAAGGACTTCAAACTCTACTCCCACT
AATAGCTTTTTGATGACTTCTAGCAAGCCTCGCTAACCTCGCCTTACCCCCCACTATTAA
CCTACTGGGAGAACTCTCTGTGCTAGTAACCACGTTCTCCTGATCAAATATCACTCTCCT
ACTTACAGGACTCAACATACTAGTCACAGCCCTATACTCCCTCTACATATTTACCACAAC
ACAATGGGGCTCACTCACCCACCACATTAACAACATAAAACCCTCATTCACACGAGAAAA
CACCCTCATGTTCATACACCTATCCCCCATTCTCCTCCTATCCCTCAACCCCGACATCAT".Replace("\n","").Replace("\r","").ToUpper();

        [TestMethod]
        public void TestVariantCalls()
        {
            var s1 = "TCTCATAATCG--CCCACGGGCTTACATCCTCATTACT-A-TTCTGCCTAGCAAACTCAAAAAAAAATTTAAAAAACT-T-".Replace("-","");
            var s2 = "TCTCATAATCGCCCCCACGGGCTTACAT--TCATTACT-G-TTCTGCCTAGCAAACTC--AAAAAAATTTAAAAAACT-A-".Replace("-","");
            var seq1 = new Sequence(DnaAlphabet.Instance, s1);
            var seq2 = new Sequence(DnaAlphabet.Instance, s2);

            var r = new Reference(seq1);
            var aln = r.AlignSequence(seq2);
            var variants = VariantCaller.VariantCaller.CallVariants(aln[0], seq1);
            Assert.AreEqual(variants.Count, 4); // Last SNP isn't called.
            
            var firstIndel = variants[0] as IndelVariant;
            Assert.AreEqual(10,firstIndel.StartPosition);
            Assert.AreEqual<IndelType>(IndelType.Insertion,firstIndel.InsertionOrDeletion);
            Assert.AreEqual("CC", firstIndel.InsertedOrDeletedBases);
            Assert.AreEqual(2, firstIndel.Length);

            var indel2 = variants[1] as IndelVariant;
            Assert.AreEqual<int>(25, indel2.StartPosition);
            Assert.AreEqual<string>("CC", indel2.InsertedOrDeletedBases);
            Assert.AreEqual<IndelType>(IndelType.Deletion, indel2.InsertionOrDeletion);
            
            var snp = variants[2] as SNPVariant;
            Assert.AreEqual(snp.StartPosition, 36);
            Assert.AreEqual<char>(snp.AltBP, 'G');
            Assert.AreEqual<int>(snp.Length, 1);

            var indel3 = variants[3] as IndelVariant;
            Assert.AreEqual<int>(53, indel3.StartPosition);


        }

        [TestMethod]
        public void TestMatch()
        {
            var r = new Sequence (DnaAlphabet.Instance, "CCCGGGGATCCTCTAGAATGCATCAGTAGAGTACGATGCTACAGCTGTGACTGTGCGCACTGCTGAGTCTGTCACTCATGTATCTGCTACGTCGTCTCACGCATACTAGACGACATGAGCGCTACGTCGAGCGTAGCAGAGATAGTATGACTCTCAGTCATATACACACATCGTGACGATGCAGAGCGATCTATCGCGCTCGCATATAGTGTGATCAAGCTTGCTGAGGACTAGTAGCTTC");
            var q = new Sequence (DnaAlphabet.Instance,    "TGGGATCCTCTAGAATGCATCAGTAGAGTACGATGCTACAGCTGTGACTGTGCGCACTGCTGAGTCTGTCAC");

            var r2 = new Reference (r);
            var a2 = r2.AlignSequence (q);
            var i = a2 [0].FindQueryPositionCorrespondingtoReferencePosition (5);
            Assert.AreEqual<int>(2, i);
        }


    }
}
