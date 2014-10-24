using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using NUnit.Framework;
using PacBio.Align;
using PacBio.HDF;
using PacBio.IO;
using PacBio.IO.Fasta;
using PacBio.Utils;

using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Single;

namespace PacBio.Consensus.Test
{
    [TestFixture]
    public class ClusteringTest
    {
        [Test]
        public void TestNativeProvider()
        {
            bool hasNativeProvider;

            // Set up the MKL Linear algebra provider, if available.
            // The MKL dll must be in the build directory for this to work.
            // using MKL makes the coarse clustering _much_ faster.
            try
            {
                Control.LinearAlgebraProvider =
                    new MathNet.Numerics.Providers.LinearAlgebra.OpenBlas.OpenBlasLinearAlgebraProvider();

                var r = new System.Random();

                // Do a test matrix multiply
                var size = 10;
                Matrix<float> m = DenseMatrix.Create(size, size, (i, j) => (float) r.NextDouble());
                var mNew = m.Multiply(m);
                mNew.Multiply(m);

                hasNativeProvider = true;
            }
            catch (Exception)
            {
                hasNativeProvider = false;
            }

            Assert.True(hasNativeProvider);
        }

        [Test]
        public void TestBinomialApprox()
        {

            var r = CoarseClustering.OccupancyBounds(5, 100, 0.95);

            var lowerBound = r.Item1;
            var upperBound = r.Item2;

            Assert.Less(lowerBound, 0.04);
            Assert.Greater(upperBound, 0.06);
        }

        [Test, Ignore]
        public void RunClustering()
        {
            Console.WriteLine("Loading bas data...");

            BasCollection bas;
            if (Environment.OSVersion.Platform != PlatformID.Win32NT)
            {
                var fn = @"/mnt/data3/vol53/2450440/0003/Analysis_Results/m130206_084201_42161_c100473711270000001823071506131392_s1_p0.bas.h5";
                bas = new BasCollection(new[] { fn });
            }
            else
            {
                var fn = @"\\usmp-data3\vol53\2450440\0003\Analysis_Results\m130206_084201_42161_c100473711270000001823071506131392_s1_p0.bas.h5";
                //var fn = @"C:\Users\pmarks\Desktop\m130206_084201_42161_c100473711270000001823071506131392_s1_p0.bas.h5";
                bas = new BasCollection(new[] { fn });
            }


            var n = 500;
            var reads = bas.Subreads.Where(sr => sr.ReadScore > 0.78 && sr.Region.Length > 2600).Take(n).ToArray();


            Console.WriteLine("Coarse Clustering reads...");
            var clusters = CoarseClustering.ClusterReads(reads, 1.6f);

            for (int i = 0; i < clusters.Length; i++)
            {
                SubreadsToFasta(clusters[i], String.Format("clusters-{0}.fasta", i));
            }
            /*
            Console.WriteLine("Fine clustering reads...");

            var clusterConsensus =
                clusters.Where(readSet => readSet.Length > 20).
                ParSelect(cluster => FineClustering.FineClusteringAndConsensus(0, cluster)).
                SelectMany(v => v);
            

            var fastaOut = @"consensus-all-new.fasta";
            File.Delete(fastaOut);

            var fr = new SimpleFASTAWriter(fastaOut);
            clusterConsensus.Apply((i, str) => fr.WriteEntry(str.FastaName, str.Sequence));
            fr.Close();
             */
        }

        [Test, Ignore]
        public void TooManyAllocated()
        {
            Console.WriteLine("Loading bas data...");

            BasCollection bas;
            if (Environment.OSVersion.Platform != PlatformID.Win32NT)
            {
                var fn = @"/mnt/data3/vol53/2450440/0003/Analysis_Results/m130206_084201_42161_c100473711270000001823071506131392_s1_p0.bas.h5";
                bas = new BasCollection(new[] { fn });
            }
            else
            {
                var fn = @"C:\Users\pmarks\Desktop\m130206_084201_42161_c100473711270000001823071506131392_s1_p0.bas.h5";
                bas = new BasCollection(new[] { fn });
            }
            
            // Example 1
            //var example = bas.Subreads.First(sr => sr.SubreadId == "m130206_084201_42161_c100473711270000001823071506131392_s1_p0/728/0_5723");
            //var reg = new AlignedSequenceReg(260, 5619, 0, 4771, Strand.Reverse);

            // Example 2
            var example = bas.Subreads.First(sr => sr.SubreadId == "m130206_084201_42161_c100473711270000001823071506131392_s1_p0/1127/0_5675");
            var reg = new AlignedSequenceReg(83, 5428, 0, 4771, Strand.Forward);
            

            var consensus = SimpleFASTAReader.ReadEntries(@"Test/tooManyAllocatedTestConsensus.fasta").First().GetSequence();

            var tpl = new TrialTemplate()
                {
                    Sequence = consensus,
                    StartAdapterBases = 0,
                    EndAdapterBases = 0
                };


            using (var qvTable = ParameterLoading.DefaultQuiver)
            {
                var s = new MultiReadMutationScorer(new[] { new Tuple<AlignedSequenceReg, IZmwBases>(reg, example.Bases) },
                                                    tpl, qvTable);

                Assert.Less(s.AllocatedEntries[0], 1000000);
            }
        }


        [Test, Explicit]
        public void MatrixMultiplySpeed()
        {
            var sizes = new[] { 100, 200, 400, 800, 1600, 3200 };
            var tf = new [] { true, false };
            int n = 1;

            //MathNet.Numerics.Control.ParallelizeOrder = 512;

            foreach (var sz in sizes)
            {
                foreach (var t in tf)
                {
                    var time = DoMatrixMultiply(sz,t, n);
                    Console.WriteLine("Size: {0}, DisableParallel: {1}, n: {2} --- {3} ms", sz, t, n, time);

                    //time = DoMatrixMultiplyStupid(sz,t, n);
                    //Console.WriteLine("Stupid -- Size: {0}, DisableParallel: {1}, n: {2} --- {3} ms", sz, t, n, time);
                }
            }


        }

        public long DoMatrixMultiply(int size, bool disablePar, int n)
        {
            // Just let it run parallel - seems to have a massive slowdown with parallelization
            //Control.LinearAlgebraProvider = new MathNet.Numerics.Algorithms.LinearAlgebra.Mkl.MklLinearAlgebraProvider();
           
            if (disablePar)
            {
                MathNet.Numerics.Control.UseSingleThread();
            }

            var r = new System.Random();

            Matrix<float> m = DenseMatrix.Create(size, size, (i, j) => (float) r.NextDouble());

            var sw = new System.Diagnostics.Stopwatch();
            sw.Start();

            for (int i = 0; i < n; i++)
            {
                // Expansion
                var mNew = m.Multiply(m);
                m = mNew;
            }

            sw.Stop();

            return sw.ElapsedMilliseconds;
        }

        public static double[,] StupidMultiply(float[,] m1, float[,] m2)                  // Stupid matrix multiplication
        {
            if (m1.GetLength(1) != m2.GetLength(0)) throw new Exception("Wrong dimensions of matrix!");

            double[,] result = new double[m1.GetLength(0), m2.GetLength(1)];
            //Matrix result = ZeroMatrix(m1.rows, m2.cols);

            var rRows = result.GetLength(0);
            var rCols = result.GetLength(1);
            var m1Cols = m1.GetLength(1);

            for (int i = 0; i < rRows; i++)
                for (int j = 0; j < rCols; j++)
                    for (int k = 0; k < m1Cols; k++)
                        result[i, j] += m1[i, k] * m2[k, j];
            return result;
        }


        public void SubreadsToFasta(IEnumerable<Subread> sr, string fn)
        {
            var fw = new SimpleFASTAWriter(fn);

            foreach (var s in sr)
            {
                fw.WriteEntry(s.SubreadId, s.FwdSequence);
            }

            fw.Close();
        }
    }
}
