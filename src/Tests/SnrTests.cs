using System;
using NUnit.Framework;
using ConsensusCore;

namespace Tests
{
    [TestFixture ()]
    public class SnrTests
    {
       
        [Test ()]
        public void testSnrCalibration ()
        {                
            SNR settings = new SNR (6.0, 6.0, 6.0, 6.0);

            // Compare to values calculated in R.
            var result = ContextParameterProvider.GetTransitionParameters ("NA", settings);
            double eps = 0.001;
            Assert.That (Math.Abs(result.Branch - -3.2866373) < eps);
            Assert.That (Math.Abs (-2.7102493 - result.Deletion) < eps);
            Assert.That (Math.Abs (-0.1375533 - result.Match) < eps);
            Assert.That (Math.Abs (-3.7044989 - result.Stick) < eps);

            //Now let's check to make sure the right channel is used.
            settings = new SNR (0.0, 8.0, 0.0, 0.0);
            result = ContextParameterProvider.GetTransitionParameters ("NC", settings);
            Assert.That (Math.Abs(result.Branch - -3.63717488) < eps);
            Assert.That (Math.Abs (-3.55678567 - result.Deletion) < eps);
            Assert.That (Math.Abs (-0.09390532 - result.Match) < eps);
            Assert.That (Math.Abs (-3.35888384 - result.Stick) < eps);
        }
    }
}

