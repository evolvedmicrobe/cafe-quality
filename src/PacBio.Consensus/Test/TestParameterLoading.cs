using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using NUnit.Framework;

namespace PacBio.Consensus.Test
{
    /// <summary>
    /// Test that we load Quiver parameter sets from the ini file correctly
    /// </summary>
    [TestFixture]
    public class TestParameterLoading
    {
        [Test]
        public void TestP4C2()
        {
            var fn = ParameterLoading.SelectParameterFile(null, "QuiverParameters.ini");
            try
            {
                using (var t = ParameterLoading.LoadParametersFromFile(fn, "P4-C2", "AllQVsMergingByChannelModel"))
                {
                    // If we get this far then we loaded the P4-C2.AllQVs.MergingByChannelModel without incident
                    Console.WriteLine("Successfully loaded 'P4-C2.AllQVsMergingByChannelModel'");
                }
            }
            catch
            {
                Assert.Fail("Failed to load 'P4-C2.AllQVsMergingByChannelModel'");
            }
        }
        
        [Test]
        public void TestUnknown()
        {
            var fn = ParameterLoading.SelectParameterFile(null, "QuiverParameters.ini");
            try
            {
                using (var t = ParameterLoading.LoadParametersFromFile(fn, "unknown"))
                {
                    Console.WriteLine("Successfully loaded 'unknown' model");
                }
            }
            catch
            {
                Assert.Fail("Failed to load 'unknown' model");
            }
        }
    }


    [TestFixture]
    public class TestChemistryMapping
    {
        [Test]
        public void TryChemistryMapping()
        {
            var m = ChemistryMapping.GetMappingForMovie("Test/chemistry_mapping.xml", "m130411_210618_42207_c100476302550000001823069506131340_s1_p0");
            Assert.AreEqual(m, "C2");
        }

    }


}
