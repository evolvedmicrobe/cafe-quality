using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Runtime.Serialization.Formatters.Binary;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Single;
using NUnit.Framework;
using PacBio.IO;
using PacBio.FSharp.Utils;

namespace PacBio.Consensus.Test
{
    /// <summary>
    /// Tests of fine clustering / phasing code
    /// </summary>
    [TestFixture]
    public class Phasing
    {
        [TestFixtureSetUp]
        public void Init()
        {
            try
            {
                Control.LinearAlgebraProvider = new MathNet.Numerics.Providers.LinearAlgebra.OpenBlas.OpenBlasLinearAlgebraProvider();
                var r = new System.Random();

                // Do a test matrix multiply
                const int size = 10;
                Matrix<float> m = DenseMatrix.Create(size, size, (i, j) => (float)r.NextDouble());
                var mNew = m.Multiply(m);
                mNew.Multiply(m);
            }
            catch
            {
                // Fall back to the managed provider
                Control.LinearAlgebraProvider = new MathNet.Numerics.Providers.LinearAlgebra.ManagedLinearAlgebraProvider();
            }
        }

        private string ConvertPath(string linuxPath)
        {
            if (Environment.OSVersion.Platform == PlatformID.Win32NT)
            {
                return linuxPath.Replace("/", "\\").Replace(@"\mnt\data3\", @"\\usmp-data3\");
            }
            else
            {
                return linuxPath;
            }
        }

        public static string NormalizeReadId(string fastaHeader)
        {
            var parts = fastaHeader.Split('/');
            var movie = parts[0];
            var holeNumber = parts[1];
            return movie + "/" + holeNumber;
        }

        /// <summary>
        /// Get read whitelists from txt files
        /// </summary>
        IEnumerable<string> GetIdsTxt(string txtFile)
        {
            var lines = System.IO.File.ReadAllLines(txtFile);
            return lines.Select(NormalizeReadId);
        }

        /// <summary>
        /// Test phasing / consensus on a canned ANRI example.
        /// </summary>
        [Test]
        public void MultiPhasingTest()
        {
            const int N = 4;

            // This test works on bas.h5 files on LIMS. 
            // These may disappear at some point. The data is marked Gold so this shouldn't happen.
            // This test can be recreated by making whitelists for different haplotypes.
            var fofn = @"Test/multi-phasing-test.fofn";
            var bas = BasCollection.FromFofn(fofn);

            // Each AN*.txt file is a whitelist of reads from one haplotype
            // You can do phasing of a different number of haplotypes by combining data from more whitelists
            var readIdFiles = Directory.EnumerateFiles(@"Test/", "AN*.txt").Take(N);
            var reads = readIdFiles.SelectMany(f =>
            {
                var ids = GetIdsTxt(f);
                var sr =
                    ids.Select(id => bas.GetSubreadsForZmw(id).OrderByDescending(r => r.Region.Length).First())
                        .Where(r => r.Region.Length > 2600 && r.Region.AdapterHitAfter && r.Region.AdapterHitAfter);

                return sr.Take(40);
            }).ToArray();

            using (var qvTable = ParameterLoading.DefaultQuiver)
            {
                var results = FineClustering.FineClusteringAndConsensus(0, reads, qvTable, "0", 300);

                // Expect 2 phases to come out
                Assert.GreaterOrEqual(results.Length, N);
            }
        }

        [Test, Explicit]
        public void CreatePhasingTestCase()
        {
            Tuple<float[][], int[]> phasingData = new Tuple<float[][], int[]>(null, null);

            var binaryFormatter = new BinaryFormatter();

            // FIXME -- if you want to create a new canned haplotyping dataset
            // you need to load some reads as in MultiPhasingTest(), then 
            // hook into the FineClustering.FineClusteringAndConsensus code,
            // run POA, then generate the Quiver score matrix as in MultiTemplateConsensus.MultiSplit()
            // then write out the float[][] mutScoreVectors and the int[] positions using the below code 

            using (
                var serializationStream = new FileStream(@"Test/new-canned-phasing-test.bin", FileMode.Create,
                    FileAccess.ReadWrite, FileShare.None))
            {
                // A canned set of mutation scores for a synthetic mixture of 6 HLA-A loci
                //Each loci has 40 reads, so if the phasing were perfect we would get 40 reads in each bucket: 
                // [n .. n+40] where n = 0..5

                binaryFormatter.Serialize(serializationStream, phasingData.Item1);
                binaryFormatter.Serialize(serializationStream, phasingData.Item2);
            }
        }


        [Test]
        public void CannedMultiPhaseTest()
        {
            float[][] quiverScoreVector;
            int[] mutationPositions;

            var binaryFormatter = new BinaryFormatter();
            using (
                var serializationStream = new FileStream(@"Test/canned-phasing-test-6haplotypes-240reads.bin", FileMode.Open,
                    FileAccess.Read,
                    FileShare.None))
            {
                // A canned set of mutation scores for a synthetic mixture of 6 HLA-A loci
                //Each loci has 40 reads, so if the phasing were perfect we would get 40 reads in each bucket: 
                // [n .. n+40] where n = 0..5
                quiverScoreVector = (float[][]) binaryFormatter.Deserialize(serializationStream);
                mutationPositions = (int[]) binaryFormatter.Deserialize(serializationStream);
            }
            var splitResult = MultiPhaseSplitter.splitHaplotypes(quiverScoreVector, mutationPositions, 6);
            var splitScore = splitResult.Item1;
            var mutsToUse = splitResult.Item2;
            var readCounts = splitResult.Item3;
            var readFractions = splitResult.Item4;
            var readPosteriors = splitResult.Item5; 

            // This test set has known good haplotying score of at least 15000
            Assert.Greater(splitScore, 15000.0);

            // This test set was synthetically created to have 6 primary haplotypes (there may be more present 
            // due to chimeras)
            Assert.GreaterOrEqual(splitResult.Item2.Length, 6);
        }
    }
}
