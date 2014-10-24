using System;
using NUnit.Framework;

namespace PacBio.Consensus.Test
{
    [TestFixture]
    public class ConsensusCoreTests
    {
        [Test]
        public void SumProductCombiner()
        {
            var worst = 0.0;

            var r = new Random();

            for (int i = 0; i < 100000; i++)
            {
                var x = (float) (20.0 * r.NextDouble() - 10.0);
                var y = (float) (20.0 * r.NextDouble() - 10.0);

                var max = Math.Max(x, y);
                var diff = Math.Min(x, y) - max;
                var correct = max + Math.Log(1 + Math.Exp(diff));

                var approx = ConsensusCore.SumProductCombiner.Combine(x, y);


                var absErr = Math.Abs(correct - approx);
                worst = Math.Max(worst, absErr);
            }

            Assert.Less(worst, 0.01);
        }
    }
}
