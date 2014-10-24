using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;

namespace PacBio.Consensus.Test
{
    [TestFixture]
    public class MutationSelectionTest
    {
        private static Func<int, int, MutationScore> dm =
            (p, s) => new MutationScore {Mutation = new Mutation {TemplatePosition = p}, Score = s};

        private List<MutationScore> muts1 = new List<MutationScore> { dm(1, 6), dm(3, 15), dm(5, 6), dm(8, 10) };
        private List<MutationScore> muts2 = new List<MutationScore> { dm(1, 4), dm(2, 25), dm(5, 6), dm(9, 10) };
        private List<MutationScore> muts3 = new List<MutationScore> { dm(1, 4), dm(2, 5), dm(3, 6), dm(9, 10), dm(20,50) };


        [Test]
        public void TestMutationSelection1()
        {
            var bm = SpacedSelector.BestMutations(muts1, 3);

            Assert.That(bm.Select(m => m.Score).Sum() == 25);
            Assert.That(bm.Count == 2);
        }

        [Test]
        public void TestMutationSelection2()
        {
            var bm = SpacedSelector.BestMutations(muts2, 3);

            Assert.That(bm.Select(m => m.Score).Sum() == 35);
            Assert.That(bm.Count == 2);
        }

        [Test]
        public void TestMutationSelection3()
        {
            var bm = SpacedSelector.BestMutations(muts3, 3);

            Assert.That(bm.Select(m => m.Score).Sum() == 66);
            Assert.That(bm.Count == 3);
        }


    }
}
