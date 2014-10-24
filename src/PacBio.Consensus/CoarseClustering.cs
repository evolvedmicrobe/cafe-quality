using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Single;
using PacBio.IO;
using PacBio.Utils;
using QuickGraph;
using QuickGraph.Algorithms.MinimumSpanningTree;
using QuickGraph.Algorithms.Search;
using QuickGraph.Graphviz;
using QuickGraph.Graphviz.Dot;
using PacBio.Align;
using PacBio.FSharp.Utils;

namespace PacBio.Consensus
{
    public class PairwiseFragments
    {
        /// <summary>
        /// Summarize the pairwise fragment alignment between two subreads.
        /// </summary>

        public static PairwiseFragments AlignAndCreate(SparseAligner al, Subread read1, Subread read2)
        {
            var s1 = read1.FwdSequence;
            var s2 = read2.FwdSequence;
            var s2Rev = DNA.ReverseComplement(s2);

            var al1 = al.SparseAlignStrings(s1, s2);
            var al2 = al.SparseAlignStrings(s1, s2Rev);

            if (al1.Count > al2.Count)
            {
                return new PairwiseFragments(al1, read1, read2);
            }
            else
            {
                return new PairwiseFragments(al2, read1, read2);
            }
        }

        public static PairwiseFragments FromAligment(List<Fragment> fragments, Subread r1, Subread r2)
        {
            return new PairwiseFragments(fragments, r1, r2);
        }

        /// <summary>
        /// A measure of the 'accuracy' of the alignment, based on the density of chained k-mer matches
        /// </summary>
        public float Accuracy { get; private set; }

        /// <summary>
        /// Average of the read length of the alignment on the two subreads
        /// </summary>
        public float Length { get; private set; }

        /// <summary>
        /// Length of the alignment in the first subread
        /// </summary>
        public int Length1 { get; private set; }

        /// <summary>
        /// Length of the alignment in the second subread
        /// </summary>
        public int Length2 { get; private set; }

        /// <summary>
        /// Number of aligned fragments
        /// </summary>
        /// <value>The count.</value>
        public int Count { get; private set; }

        /// <summary>
        /// Subread 1 in alignment
        /// </summary>
        public Subread Read1 { get; private set; }

        /// <summary>
        /// Subread 2 in alignment
        /// </summary>
        public Subread Read2 { get; private set; }

        /// <summary>
        /// Fraction of the bases in the two subreads that falls inside this alignment
        /// </summary>
        public float OverlapFraction { get; private set; }

        private PairwiseFragments(List<Fragment> fragments, Subread s1, Subread s2)
        {
            var f = fragments[0];
            var l = fragments[fragments.Count - 1];

            Count = fragments.Count;
            Length = ((l.X - f.X) + (l.Y - f.Y))/2;
            Accuracy = ((float) fragments.Count)/Length;

            // Compute overlap fraction -- fraction of bases in both reads that are within the alignment
            Length1 = l.X - f.X;
            Length2 = l.Y - f.Y;
            //OverlapFraction = ((float)(Length1 + Length2))/ ((float)(s1.Region.Length + s2.Region.Length));
            OverlapFraction = Math.Max(((float)Length1) / ((float)s1.Region.Length), ((float)Length2) / ((float)s2.Region.Length));

            Read1 = s1;
            Read2 = s2;
        }
    }
    
    /// <summary>
    /// Perform a (hopefully fast) 'coarse' clustering of a set of subreads.  
    /// Coarse clustering should reads coming from templates that differ by > ~4% into different clusters.
    /// 
    /// HACKING NOTE:  we had discussed an improved version of this code which split clusters more aggresively, 
    /// then re-merged clusters if their consensus sequences (from POA) were less than x% different where x might equal 3. This would 
    /// give a more 'explainable' clustering.  Unclear how valuable that might be.
    /// </summary>
    public class CoarseClustering
    {
        private static PacBioLogger logger = PacBioLogger.GetLogger("CoarseClustering");

        private static void Log(LogLevel level, string msg)
        {
            logger.Log(level, msg);
        }

        private static void Log(LogLevel level, string msg, params object[] args)
        {
            logger.Log(level, String.Format(msg, args));
        }


        static CoarseClustering()
        {
            // Set up the MKL Linear algebra provider, if available.
            // The MKL dll must be in the build directory for this to work.
            // using MKL makes the coarse clustering _much_ faster.
            try 
            {
                Control.LinearAlgebraProvider = new MathNet.Numerics.Providers.LinearAlgebra.OpenBlas.OpenBlasLinearAlgebraProvider();
                var r = new System.Random();

                // Do a test matrix multiply
                const int size = 10;
                Matrix<float> m = DenseMatrix.Create(size, size, (i, j) => (float) r.NextDouble());
                var mNew = m.Multiply(m);
                mNew.Multiply(m);

                Log(LogLevel.DEBUG, "Using Intel MKL Linear Algebra.");
            }
            catch (Exception e)
            {
                Log(LogLevel.WARN, "Error initializing OpenBLAS native provider. Falling back to Managed Linear Algebra. Coarse clustering may be slow for large numbers of reads.");
                Log(LogLevel.WARN, "OpenBLAS error message: '{0}'", e.Message);

                // Fall back to the managed provider
                Control.LinearAlgebraProvider = new MathNet.Numerics.Providers.LinearAlgebra.ManagedLinearAlgebraProvider();
            }
        }


        public static ushort K = 10;

        // This could be explored for better overlap graph
        public static int MaxNBest = 12;

        /// <summary>
        /// Use a suffix array of the reads to find the NBest overlapping reads for each read input read.
        /// Return the pairwise overlap matrix, which contains null values for non-overlapping reads.
        /// </summary>
        public static PairwiseFragments[,] SuffixArrayOverlapFinder(Subread[] reads)
        {
            var n = reads.Length;
            var alignments = new PairwiseFragments[n, n];

            var idx = SuffixArray.SuffixArrayIndex.MakeSubreadIndex(reads);

            var nBest = (int) Math.Min(MaxNBest, Math.Sqrt(reads.Length) / 2);

			// Process the overlaps lazily to avoid accumulating all the fragment arrays in memory at the same time
            var results = Enumerable.Range(0, n).ParSelect(
                i =>
                {
                    var hits = Overlapper.findOverlaps(K, nBest, 40, idx, i);

                    var fragments = hits.Map(
                        hit =>
                        {
                            var target = (int) hit.Item1;
                            var overlaps = hit.Item2;
                            var pwf = PairwiseFragments.FromAligment(overlaps, reads[i], reads[target]);
                            return Tuple.Create(target, pwf);
                        });

                    return Tuple.Create(i, fragments);
                });

            foreach(var r in results)
            {                
                var src = r.Item1;
                var pairwiseHits = r.Item2;

                foreach (var hit in pairwiseHits)
                {
                    var target = hit.Item1;
					var pwf = hit.Item2;

                    alignments[src, target] = pwf;
                    alignments[target, src] = pwf;
                }
            }

            return alignments;
        }
        
        /// <summary>
        /// Scores a subread on the basis of it's completeness and non-chimeric-ness
        /// </summary>
        public static float SubreadScore(Subread sr)
        {
            //var hitBefore = sr.Region.AdapterHitBefore ? 1.0f : 0.5f;
            //var hitAfter = sr.Region.AdapterHitAfter ? 1.0f : 0.5f;

            // Do a sparse alignment of the reads to get the alignment band
            var sparse = new SparseAligner(6);

            var seq = sr.FwdSequence;
            var rseq = seq.ReverseComplement();

            var rcFragments = sparse.SparseAlignStrings(seq, rseq);

            var chimeraScore = 1.0f;
            if (rcFragments.Count > 0)
            {
                var fragFirst = rcFragments[0];
                var fragLast = rcFragments[rcFragments.Count - 1];
                var alnLength = fragLast.X - fragFirst.X;

                // The more of the subread that is self-complementary the bigger the penalty
                chimeraScore = (float)Math.Pow(1.0 - ((float)alnLength / (float)seq.Length), 2.0);
            }

            var score = (float)seq.Length * chimeraScore * ((float)sr.ReadScore / 1000);

            return score;
        }


        /// <summary>
        /// Do a coarse clustering of subreads, using the Markov Clustering algorithm. Optionally write the before and after overlap graphs
        /// out in .dot format.  The subreads in each cluster are sorted using the PageRank algorithm -- this puts the most high quality and longest reads 
        /// first, which makes the POA step go more smoothly.
        /// </summary>
        public static Subread[][] ClusterReads(Subread[] reads, float inflationParameter, string graphWritePath = null, bool doClustering = true)
        {
            if (reads.Length == 0)
                return new Subread[][] { };

            // Fast all-read overlapper using suffix array
            var alignments = SuffixArrayOverlapFinder(reads);

            // Distance function
            Func<PairwiseFragments, float> distance =
                v =>
                {
                    if (v == null || v.Accuracy < 0.07)
                        return 0.0f;

                    // NOTE! if the overlap fraction metric is causing this to separate 5' from 3' chunks -- look here
                    // Brett may be interested
                    return v.Count * v.OverlapFraction;
                };

            #if false
            using (var f = new StreamWriter("pairwise-data.csv"))
            {
                f.WriteLine("idx1,idx2,count,length1,length2,overlap");
                for (int i = 0; i < alignments.GetLength(0); i++)
                {
                    for (int j = 0; j < alignments.GetLength(1); j++)
                    {
                        var pwf = alignments[i, j];
                        if (pwf != null)
                            f.WriteLine(String.Format("{0},{1},{2},{3},{4},{5}",i,j,pwf.Count,pwf.Length1,pwf.Length2,pwf.OverlapFraction));
                    }
                }
                f.Close();
            }
            #endif

            var distanceMatrix = alignments.Map(distance);

            int[][] clusters = null;

            if (doClustering)
            {
                // Do Markov clustering of the distance matrix
                var m = MarkovClustering(distanceMatrix, inflationParameter);

                // Form the graph implied by the clustering
                var clusterGraph = MatrixToGraph(m, 0.0001f);

                if (graphWritePath != null)
                {
                    var initialGraph = MatrixToGraph(distanceMatrix);

                    var targetFile = Path.Combine(graphWritePath, "initial-clustering.dot");
                    WriteUndirectedGraph(initialGraph, targetFile);

                    targetFile = Path.Combine(graphWritePath, "markov-clustering.dot");
                    WriteUndirectedGraph(clusterGraph, targetFile);
                }
                // Get the connected components of the graph -- these are the clusters
                clusters = GetComponents(clusterGraph).ToArray();
            }
            else
            {
                // If the clustering is turned off, just put everything into 1 cluster
                clusters = new[] {Enumerable.Range(0, reads.Length).ToArray()};
            }

            // Do a page rank on all the reads -- use the global
            // rank to order reads in each cluster -- makes the POA process go much more smoothly
            var initialScores = reads.ParSelect(SubreadScore);
            var pageRankScores = PageRank(distanceMatrix, initialScores);

            var reorderedClusters = clusters.Map(cluster =>
            {
                var orderedCluster = cluster.OrderByDescending(readIdx => pageRankScores[readIdx]).ToArray();
                var clusterReads = orderedCluster.Map(readIdx => reads[readIdx]);  
                return clusterReads;
            });

            return reorderedClusters;
        }

        /// <summary>
        /// Convert a weight matrix to an undirected graph
        /// </summary>
        public static UndirectedGraph<int, Edge<int>> MatrixToGraph(float[,] m, float threshold = 0.0001f, Func<int, bool> nodePred = null)
        {
            var n = m.GetLength(0);
            nodePred = nodePred ?? ((v) => true);

            var graph = new UndirectedGraph<int, Edge<int>>();
            for (int i = 0; i < n; i++)
            {
                if (nodePred(i))
                    graph.AddVertex(i);
            }

            for (int i = 0; i < n; i++)
            {
                for (int j = i + 1; j < n; j++)
                {
                    if (nodePred(i) && nodePred(j) && (m[i, j] > threshold || m[j, i] > threshold))
                    {
                        var e = new Edge<int>(i, j);
                        graph.AddEdge(e);
                    }
                }
            }

            return graph;
        }

        /// <summary>
        /// Write an undirected graph to a dot file
        /// </summary>
        public static void WriteUndirectedGraph(UndirectedGraph<int, Edge<int>> graph, string filename)
        {
            var graphviz = new GraphvizAlgorithm<int, Edge<int>>(graph);
            graphviz.FormatVertex += Format;
            graphviz.FormatEdge += (o, e) =>
            {
                var edge = e.Edge;
                e.EdgeFormatter.Dir = GraphvizEdgeDirection.None;
            };

            // render
            string output = graphviz.Generate();
            System.IO.File.WriteAllText(filename, output);
        }
        
        public static void Format(Object sender, FormatVertexEventArgs<int> e)
        {
            e.VertexFormatter.Label = e.Vertex.ToString();
        }

        /// <summary>
        /// Pick out the largest components of the clustered graph
        /// </summary>
        public static TA[][] GetComponents<TA, TB>(UndirectedGraph<TA, TB> graph) where TB : IEdge<TA>
        {
            var cc = new QuickGraph.Algorithms.ConnectedComponents.ConnectedComponentsAlgorithm<TA, TB>(graph);
            cc.Compute();

            var groups = cc.Components.GroupBy(kvp => kvp.Value, kvp => kvp.Key);

            var longestComponents = groups.Select(g => g.ToArray()).OrderByDescending(c => c.Length).Where(g => g.Length > 3).ToArray();
            return longestComponents;
        }


        /// <summary>
        /// Pick out the largest components of the clustered graph
        /// </summary>
        public static QuickGraph.Algorithms.ConnectedComponents.ConnectedComponentsAlgorithm<TA, TB> GetComponents2<TA, TB>(UndirectedGraph<TA, TB> graph) where TB : IEdge<TA>
        {
            var cc = new QuickGraph.Algorithms.ConnectedComponents.ConnectedComponentsAlgorithm<TA, TB>(graph);
            cc.Compute();

            return cc;
        }


        /// <summary>
        /// Perform Markov clustering on an undirected graph.  Return the clustered adjacency matrix.
        /// See http://www.micans.org/mcl/ for more details on algorithm.
        /// </summary>
        public static float[,] MarkovClustering(float[,] similarityMatrix, double inflation = 1.6)
        {
            Matrix<float> m = DenseMatrix.OfArray(similarityMatrix);

            m.SetDiagonal(DenseVector.Create(m.RowCount, i => 1.0f));
            m = m.NormalizeColumns(1);

            for (int i = 0; i < 16; i++)
            {
                // Expansion
                var mNew = m.Multiply(m);

                // Inflation
                mNew.MapInplace(v => (float) Math.Pow(v, inflation));
                mNew = mNew.NormalizeColumns(1);
                 
                // See how much we have changed
                var diffMat = mNew - m;

                // Comment this out because it's not supported by our hacked MKL
                //var diff = (mNew - m).InfinityNorm();
                var diff = InfinityNorm(diffMat);

                m = mNew;

                if (diff < 0.00001)
                    break;
            }

            return m.ToArray();
        }

        /// <summary>
        /// Find the maximal value in a matrix -- just a hack until InfinityNorm is fixed
        /// </summary>
        static float InfinityNorm(Matrix<float> m)
        {
            var rows = m.RowCount;
            var cols = m.ColumnCount;

            var v = float.NegativeInfinity;

            for (int r = 0; r < rows; r++)
            {
                for (int c = 0; c < cols; c++)
                {
                    v = Math.Max(v, Math.Abs(m.At(r, c)));
                }
            }

            return v;
        }


        /// <summary>
        /// Rank the reads in the cluster using a PageRank style ranking.
        /// Each read starts with a score determined by it's Read Score and whether it has
        /// adapters on each end.  Then the 'endorsements' from overlapping reads contribute to the ranking.
        /// We will use this ranking to put reads into the POA in a good order.
        /// </summary>
        public static float[] PageRank(float[,] endorsementMatrix, IEnumerable<float> initialDistribution, float damping = 0.9f)
        {
            // Endorsement matrix
            Matrix<float> m = DenseMatrix.OfArray(endorsementMatrix);
            m = m.NormalizeColumns(1);

            // Initial quality of reads
            Vector<float> initialVect = DenseVector.OfEnumerable(initialDistribution);
            initialVect = initialVect.Normalize(1);

            var rankVect = initialVect.Clone();

            for (int i = 0; i < 5; i++)
            {
                rankVect = m.Multiply(rankVect);
                rankVect = rankVect * damping + (1 - damping) * initialVect;
            }

            return rankVect.ToArray();
        }



 
        /// <summary>
        /// Compute the confidence interval of the probability of a cluster. Used for pruning clusters
        /// </summary>
        public static Tuple<double, double> OccupancyBounds(int nObs, int nTotal, double confidence)
        {
            // See http://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval

            var norm = new MathNet.Numerics.Distributions.Normal(0, 1.0);
            var z = norm.InverseCumulativeDistribution(1 - (1 - confidence)/2);

            var pHat = ((double) nObs)/nTotal;

            var lb = (pHat + z*z/(2*nTotal) - z*Math.Sqrt((pHat*(1 - pHat) + z*z/(4.0*nTotal))/nTotal))/
                     (1.0 + z*z/nTotal);

            var ub = (pHat + z*z/(2*nTotal) + z*Math.Sqrt((pHat*(1 - pHat) + z*z/(4.0*nTotal))/nTotal))/
                     (1.0 + z*z/nTotal);

            return new Tuple<double, double>(lb, ub);
        }
    }
}
