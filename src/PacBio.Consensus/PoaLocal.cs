using System;
using System.Collections;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using PacBio.Align;
using PacBio.IO;
using PacBio.Utils;
using QuickGraph;
using QuickGraph.Algorithms;
//using QuickGraph.Graphviz;
using QuickGraph.Graphviz;
using QuickGraph.Graphviz.Dot;
using PartialOrderPath = PacBio.Consensus.PartialOrderPath<PacBio.Consensus.AVertex>;

namespace PacBio.Consensus
{
    
    /// <summary>
    /// An implementation of Partial Order Alignments, as described in 'Mulitple sequence alignment 
    /// using partial order graphs', Lee, Grasso and Sharlow, Bioinformatics, vol. 18, no. 3, 2002, pg 452.
    /// 
    /// Implements the alignment between a partial order graph and a sequence. The alignment is used to update 
    /// the partial order graph with the new sequence.  A simple max-score recursion is used to compute a 
    /// consensus sequence from the partial-order graph.
    /// 
    /// This code is intended to be a very fast initial guess at the consensus of a group of reads.
    /// </summary>
    public class PoaLocal
    {
        /// <summary>
        /// Manage arrays used during POA alignment 
        /// </summary>
        private class BufferProvider
        {
            private List<Tuple<int[], PartialOrderPath[]>> arraysInUse;

            internal BufferProvider()
            {
                arraysInUse = new List<Tuple<int[], PartialOrderPath[]>>();
            }

            /// <summary>
            /// Allocate a score and path array.  These may have old values in them -- caller 
            /// </summary>
            /// <param name="length"></param>
            /// <returns></returns>
            public Tuple<int[], PartialOrderPath[]> GetArrays(int length)
            {
                int roundSize = MathUtils.PowUp(length + 1, 2);

                Tuple<int[], PartialOrderPath[]> t;
                var queue = GetQueue(roundSize);

                if (!queue.TryDequeue(out t))
                {
                    t = new Tuple<int[], PartialOrderPath[]>(new int[roundSize], new PartialOrderPath[roundSize]);
                }
                else
                {
                    // The arrays have already been used -- clean them up
                    var score = t.Item1;
                    var path = t.Item2;
                    
                    for (int i = 0; i < score.Length; i++)
                    {
                        score[i] = 0;
                        path[i] = new PartialOrderPath();
                    }
                }

                arraysInUse.Add(t);
                return t;
            }

            /// <summary>
            /// Return all arrays in use to the pool
            /// </summary>
            public void ReturnAllArrays()
            {
                foreach (var r in arraysInUse)
                {
                    var sz = r.Item1.Length;
                    var queue = GetQueue(sz);
                    queue.Enqueue(r);
                }

                arraysInUse.Clear();
            }
        }

        /// <summary>
        /// Maintain a queue of ununsed buffers for each power of 2 size.  GetQueue get the queue of the correct size
        /// </summary>
        private static ConcurrentQueue<Tuple<int[], PartialOrderPath[]>> GetQueue(int size)
        {
            return BufferSet.GetOrAdd(size, sz => new ConcurrentQueue<Tuple<int[], PartialOrderPath[]>>());
        }

        /// <summary>
        /// Shared dictionary of queues where the not-in-use arrays are kept.
        /// </summary>
        private static ConcurrentDictionary<int, ConcurrentQueue<Tuple<int[], PartialOrderPath[]>>> BufferSet =
            new ConcurrentDictionary<int, ConcurrentQueue<Tuple<int[], PartialOrderPath[]>>>();

        public static string[] FreeReport()
        {
            return BufferSet.Select(b => String.Format("Size: {0}, Free:{1}", b.Key, b.Value.Count)).ToArray();
        }
        
        /// <summary>
        /// POA Graph
        /// </summary>
        public PoaGraph<AVertex, Edge<AVertex>> Graph { get; private set; }
        private BufferProvider buffers;

        /// <summary>
        /// Length of reads in POA
        /// </summary>
        public List<int> ReadLengths;

        /// <summary>
        /// Total number of reads in alignment
        /// </summary>
        public int NumReads
        {
            get { return ReadLengths.Count; }
        }

        /// <summary>
        /// Construct the basic partial order graph of a single read
        /// </summary>
        public PoaLocal()
        {
            buffers = new BufferProvider();
            Graph = new PoaGraph<AVertex, Edge<AVertex>>(false);
            ReadLengths = new List<int>();
        }

        /// <summary>
        /// Setup node for the first read in the POA
        /// </summary>
        private void AddFirstSeq(string seq1, bool isRollIn = false)
        {
            var startVertex = new AVertex(seq1[0], NumReads, 0);
            Graph.AddVertex(startVertex);

            for (int i = 1; i < seq1.Length; i++)
            {
                var v = new AVertex(seq1[i], NumReads, i);
                v.DoNotCall = isRollIn;

                var e = new Edge<AVertex>(startVertex, v);

                Graph.AddVertex(v);
                Graph.AddEdge(e);

                startVertex = v;
            }

            ReadLengths.Add(seq1.Length);
            graphInitialized = true;
        }

        // Alignment scoring weights
        public short[] Scores = new short[] { -4, -4, -6, 3, -3 };

        /// <summary>
        /// Initial read has been added to the graph
        /// </summary>
        private bool graphInitialized = false;
        
        /// <summary>
        /// Tack on an extra column to the best columnm - this is the 'n+1' column of a SW alignment.
        /// We will start the traceback from this new vertex
        /// </summary>
        private AVertex MakeTerminal(AVertex terminus, Dictionary<AVertex, PoaColumn> dictionary, string read, bool freeInsertions)
        {
            var tt = new Tuple<AVertex, PoaColumn>(terminus, dictionary[terminus]);
            var inList = new []{ tt };

            // set up a dummy vertex that is the terminus of the alignment
            var dummyVertex = new AVertex('N');
            int d1, d2;

            var terminalColumn = AlignCol(read, inList, dummyVertex, new Range(0, read.Length), freeInsertions, out d1, out d2, Scores[3]);
            dictionary[dummyVertex] = terminalColumn;
            return dummyVertex;
        }

        /// <summary>
        /// Add successive reads to the POA
        /// </summary>
        public void AddReads(string[] reads)
        {
            // the backbone performs sparse 'fragment-wise' pairwise alignment 
            // which is used to band the POA alignment
            //var backbone = new AlignmentBackbone(reads);

            foreach (var read in reads)
            {
                AddReadSparse2(read);
            }
        }


        /// <summary>
        /// Add a new read to the partial order graph.  First compute the alignment between the partial order graph and the new read, then
        /// trace back this alignment and thread the new read into the graph.
        /// </summary>
        public void AddReadSparse(string read, AlignmentBackbone backbone)
        {
            // Add the first read
            if (!graphInitialized)
            {
                AddFirstSeq(read);
                return;
            }

            // A column of the alignment will be computed for each vertex in the POA
            var dictionary = new Dictionary<AVertex, PoaColumn>();

            AVertex curVertex = null;
            int bestReadPos = 0;

            try
            {
                int bestScore = int.MinValue;
                var topo = Graph.TopologicalSort;

                // rangeFunc returns the read range to align for a given index of the POA topological sort
                var rangeFunc = backbone.PrepareRanges(NumReads, Graph);

                for (int idx = 0; idx < topo.Length; idx++)
                {
                    var v = topo[idx];
                    // Get the alignment columns for predeccsor vertices
                    var inBoundScores =
                        Graph.InEdges(v).Select(
                            inV => new Tuple<AVertex, PoaColumn>(inV.Source, dictionary[inV.Source])).
                              ToArray();

                    bool isTerminalVertex = !Graph.OutEdges(v).Any();

                    // Make an alignment column for the current vertex, given the columns of an predecessor vertices.
                    int colBestScore;
                    int colBestReadPos;

                    var range = rangeFunc(idx);
                    var col = AlignCol(read, inBoundScores, v, range, false, out colBestScore, out colBestReadPos, Scores[3]);

                    // Cache the alignment col
                    dictionary[v] = col;

                    // Keep track of vertices the terminate the DAG

                    // Case 1 - read goes beyond end of graph
                    if (isTerminalVertex)
                    {
                        var dummyVertex = MakeTerminal(v, dictionary, read, true);
                        var dummyCol = dictionary[dummyVertex];

                        if (dummyCol.Score(read.Length) > bestScore)
                        {
                            curVertex = dummyVertex;
                            bestScore = dummyCol.Score(read.Length);
                            bestReadPos = read.Length;
                        }
                    }
                    else
                    {
                        // Scan for best score
                        if (colBestScore > bestScore)
                        {
                            curVertex = v;
                            bestScore = colBestScore;
                            bestReadPos = colBestReadPos;
                        }
                    }
                }

                AVertex nextV = null;
                int readIdx = bestReadPos;

                // Prevent chimeras from forming in doubly loaded ZMWs -- different fragments must overlap reasonably well to be put into the same chunk.
                // If the alignment score isn't big enough, this chunks goes in as an unconnected component, then we bail
                const int bestScoreThreshold = 30;
                if (bestScore < bestScoreThreshold)
                {
                    var startVertex = new AVertex(read[0], NumReads, 0);
                    Graph.AddVertex(startVertex);

                    for (int i = 1; i < read.Length; i++)
                    {
                        var v = new AVertex(read[i], NumReads, i);
                        var e = new Edge<AVertex>(startVertex, v);

                        Graph.AddVertex(v);
                        Graph.AddEdge(e);

                        startVertex = v;
                    }

                    ReadLengths.Add(read.Length);
                    return;
                }

                // If we terminated the alignment before a terminal vertex, we need to add the remaining read bases into the sequence
                for (int i = read.Length; i > bestReadPos; i--)
                {
                    var newV = new AVertex(read[i - 1], NumReads, i - 1);
                    Graph.AddVertex(newV);
                    if (nextV != null)
                    {
                        var e = new Edge<AVertex>(newV, nextV);
                        Graph.AddEdge(e);
                    }

                    nextV = newV;
                }


                // Now we have an alignment between the initial POA and the current read.  Trace it back and mutate the POA to add the new read.
                while (curVertex != null)
                {
                    var cvScores = dictionary[curVertex];
                    var pathSeg = cvScores.Path(readIdx);

                    var srcVertex = pathSeg.PrevVertex;

                    if (pathSeg.Mode == AlignMode.Insert)
                    {
                        // We make a new Vertex to represent the insertion
                        // It will rejoin the graph the next time we have a match
                        var newV = new AVertex(read[readIdx - 1], NumReads, readIdx - 1);
                        Graph.AddVertex(newV);

                        if (nextV != null)
                        {
                            var e = new Edge<AVertex>(newV, nextV);
                            Graph.AddEdge(e);
                        }

                        nextV = newV;

                        readIdx--;
                    }
                    else if (pathSeg.Mode == AlignMode.Match)
                    {
                        // Update the current vertex
                        srcVertex.AddSequence(NumReads, readIdx - 1);

                        if (nextV != null)
                        {
                            var e = new Edge<AVertex>(srcVertex, nextV);
                            Graph.AddEdge(e);
                        }

                        nextV = srcVertex;

                        readIdx--;
                    }
                    else if (pathSeg.Mode == AlignMode.Mismatch)
                    {
                        var newV = new AVertex(read[readIdx - 1], NumReads, readIdx - 1);
                        Graph.AddVertex(newV);

                        if (nextV != null)
                        {
                            var e = new Edge<AVertex>(newV, nextV);
                            Graph.AddEdge(e);
                        }

                        nextV = newV;
                        readIdx--;
                    }
                    else if (pathSeg.Mode == AlignMode.Delete)
                    {
                        // No new vertex, just update the current vertex
                    }
                    // Update to the new current vertex;
                    curVertex = pathSeg.PrevVertex;
                }

                // If the alignment started part way into the read, add the preceding bases of the read at the beginning
                while (readIdx > 0)
                {
                    var newV = new AVertex(read[readIdx - 1], NumReads, readIdx - 1);
                    Graph.AddVertex(newV);

                    var e = new Edge<AVertex>(newV, nextV);
                    Graph.AddEdge(e);

                    nextV = newV;
                    readIdx--;
                }

                ReadLengths.Add(read.Length);
            }
            finally
            {
                if (buffers != null)
                    buffers.ReturnAllArrays();
            }
        }



        public class AlignmentSummary : IDisposable
        {
            public int Score;
            public string Read;
            public int BestReadPos;
            public AVertex CurrentVertex;
            public Dictionary<AVertex, PoaColumn> PoaColumns;

            internal PoaLocal PoaGraph;

            public void Commit()
            {
                PoaGraph.CommitAdd(this);
            }

            public void Dispose()
            {
                if (PoaGraph.buffers != null)
                    PoaGraph.buffers.ReturnAllArrays();
            }
        }


        /// <summary>
        /// Call the consensus of the existing graph, do a sparse alignment between the consensus and the new read,
        /// then align the new read using the backbone provided by the sparse alignment.
        /// </summary>
        public void AddReadSparse2(string read)
        {
            using (var al = TryAddRead(read))
            {
                if (al != null)
                    al.Commit();
            }
        }

        /// <summary>
        /// Align the read to the current POA, using the consensus of the current POA as a backbone band to align to.
        /// Return an AlignmentSummary that can be inspected by calling code, then .Commit() to add the alignment to the POA,
        /// or  .Dispose() to forget about that alignment
        /// </summary>
        public AlignmentSummary TryAddRead(string read)
        {
            // Add the first read
            if (!graphInitialized)
            {
                AddFirstSeq(read);
                return null;
            }

            // A column of the alignment will be computed for each vertex in the POA
            var dictionary = new Dictionary<AVertex, PoaColumn>();

            AVertex curVertex = null;
            int bestReadPos = 0;
            int bestScore = int.MinValue;

            try
            {
                ComputeCoverage();

                var minCoverage = 1;

                Func<AVertex, float> score = v =>
                    {
                        if (v.DoNotCall || v.Coverage < minCoverage)
                            return Single.NegativeInfinity;

                        // penalize regions with low coverage -- coverage should only be able to drop to
                        // NumReads - 2 in a SMRTBell read.
                        return 2*v.NumSequences - 1*v.Coverage - 0.1f;
                    };

                var path = GraphAlgos.MaxPath(Graph, score).ToArray();
                var consensusSeq = new string(path.Select(v => v.Base).ToArray());

                var vertices = Graph.TopologicalSort;

                // rangeFunc returns the read range to align for a given index of the POA topological sort
                var rangeFunc = AlignmentBackbone.SimpleRange(read, consensusSeq, Graph, path);

                for (int idx = 0; idx < vertices.Length; idx++)
                {
                    var v = vertices[idx];

                    // Get the alignment columns for predeccsor vertices
                    var inBoundScores =
                        Graph.InEdges(v).Select(
                            inV => new Tuple<AVertex, PoaColumn>(inV.Source, dictionary[inV.Source])).
                              ToArray();

                    bool isTerminalVertex = !Graph.OutEdges(v).Any();

                    // Make an alignment column for the current vertex, given the columns of an predecessor vertices.
                    int colBestScore;
                    int colBestReadPos;

                    var range = rangeFunc(idx);
                    
                    // Compute a score for this vertex based on it's quality
                    short match = Scores[3];
                    var vertexReward = (short) Math.Min(match, Math.Round(match*((1.0*v.NumSequences/(Math.Max(1, v.Coverage - 1))))));
                    var col = AlignCol(read, inBoundScores, v, range, false, out colBestScore, out colBestReadPos, vertexReward);

                    // Cache the alignment col
                    dictionary[v] = col;

                    // Keep track of vertices the terminate the DAG

                    // Case 1 - read goes beyond end of graph
                    if (isTerminalVertex)
                    {
                        var dummyVertex = MakeTerminal(v, dictionary, read, true);
                        var dummyCol = dictionary[dummyVertex];

                        if (dummyCol.Score(read.Length) > bestScore)
                        {
                            curVertex = dummyVertex;
                            bestScore = dummyCol.Score(read.Length);
                            bestReadPos = read.Length;
                        }
                    }
                    else
                    {
                        // Scan for best score
                        if (colBestScore > bestScore)
                        {
                            curVertex = v;
                            bestScore = colBestScore;
                            bestReadPos = colBestReadPos;
                        }
                    }
                }

            }
            catch
            {
                buffers.ReturnAllArrays();
                throw;
            }

            return new AlignmentSummary
                {
                    BestReadPos = bestReadPos,
                    Read = read,
                    Score = bestScore,
                    CurrentVertex = curVertex,
                    PoaColumns = dictionary,
                    PoaGraph = this
                };
        }

        public void CommitAdd(AlignmentSummary al)
        {
            try
            {
                var bestReadPos = al.BestReadPos;
                var read = al.Read;
                var curVertex = al.CurrentVertex;
                var poaColumns = al.PoaColumns;

                AVertex nextV = null;
                int readIdx = bestReadPos;

                /*
                // Prevent chimeras from forming in doubly loaded ZMWs -- different fragments must overlap reasonably well to be put into the same chunk.
                // If the alignment score isn't big enough, this chunks goes in as an unconnected component, then we bail
                var bestScoreThreshold = 30;
                if (bestScore < bestScoreThreshold)
                {
                    var startVertex = new AVertex(read[0], NumReads, 0);
                    Graph.AddVertex(startVertex);

                    for (int i = 1; i < read.Length; i++)
                    {
                        var v = new AVertex(read[i], NumReads, i);
                        var e = new Edge<AVertex>(startVertex, v);

                        Graph.AddVertex(v);
                        Graph.AddEdge(e);

                        startVertex = v;
                    }

                    NumReads++;
                    return;
                }
                */


                // If we terminated the alignment before a terminal vertex, we need to add the remaining read bases into the sequence
                for (int i = read.Length; i > bestReadPos; i--)
                {
                    var newV = new AVertex(read[i - 1], NumReads, i - 1);
                    Graph.AddVertex(newV);
                    if (nextV != null)
                    {
                        var e = new Edge<AVertex>(newV, nextV);
                        Graph.AddEdge(e);
                    }

                    nextV = newV;
                }


                // Now we have an alignment between the initial POA and the current read.  Trace it back and mutate the POA to add the new read.
                while (curVertex != null)
                {
                    var cvScores = poaColumns[curVertex];
                    var pathSeg = cvScores.Path(readIdx);

                    var srcVertex = pathSeg.PrevVertex;

                    if (pathSeg.Mode == AlignMode.Insert)
                    {
                        // We make a new Vertex to represent the insertion
                        // It will rejoin the graph the next time we have a match
                        var newV = new AVertex(read[readIdx - 1], NumReads, readIdx - 1);
                        Graph.AddVertex(newV);

                        if (nextV != null)
                        {
                            var e = new Edge<AVertex>(newV, nextV);
                            Graph.AddEdge(e);
                        }

                        nextV = newV;

                        readIdx--;
                    }
                    else if (pathSeg.Mode == AlignMode.Match)
                    {
                        // Update the current vertex
                        srcVertex.AddSequence(NumReads, readIdx - 1);

                        if (nextV != null)
                        {
                            var e = new Edge<AVertex>(srcVertex, nextV);
                            Graph.AddEdge(e);
                        }

                        nextV = srcVertex;

                        readIdx--;
                    }
                    else if (pathSeg.Mode == AlignMode.Mismatch)
                    {
                        var newV = new AVertex(read[readIdx - 1], NumReads, readIdx - 1);
                        Graph.AddVertex(newV);

                        if (nextV != null)
                        {
                            var e = new Edge<AVertex>(newV, nextV);
                            Graph.AddEdge(e);
                        }

                        nextV = newV;
                        readIdx--;
                    }
                    else if (pathSeg.Mode == AlignMode.Delete)
                    {
                        // No new vertex, just update the current vertex
                    }
                    // Update to the new current vertex;
                    curVertex = pathSeg.PrevVertex;
                }
                

                // If the alignment started part way into the read, add the preceding bases of the read at the beginning
                while (readIdx > 0)
                {
                    var newV = new AVertex(read[readIdx - 1], NumReads, readIdx - 1);
                    Graph.AddVertex(newV);

                    if (nextV != null)
                    {
                        var e = new Edge<AVertex>(newV, nextV);
                        Graph.AddEdge(e);
                    }

                    nextV = newV;
                    readIdx--;
                }
            }
            finally
            {
                if (buffers != null)
                    buffers.ReturnAllArrays();
            }

            ReadLengths.Add(al.Read.Length);
        }

        /// <summary>
        /// Hold a banded column of the POA alignment.  Track the score matrix, the path matrix, and the start positiion of the band
        /// </summary>
        public struct PoaColumn
        {
            /// <summary>
            /// Start Base of the band
            /// </summary>
            public int StartBase;

            /// <summary>
            /// Alignment score of this column
            /// </summary>
            public int[] ScoreChunk;

            /// <summary>
            /// Alignment path of this column
            /// </summary>
            public PartialOrderPath[] PathChunk;
            
            /// <summary>
            /// Score in this column at read position p.
            /// </summary>
            public int Score(int p)
            {
                return ScoreChunk[p - StartBase];
            }

            /// <summary>
            /// Alignment path pointer at read position p.
            /// </summary>
            /// <param name="p"></param>
            /// <returns></returns>
            public PartialOrderPath Path(int p)
            {
                return PathChunk[p - StartBase];
            }

            public bool HasPosition(int p)
            {
                return p >= StartBase && p < StartBase + ScoreChunk.Length;
            }

        }

        /// <summary>
        /// Make one column of the POA alignment. 
        /// </summary>
        /// <param name="read">The read string</param>
        /// <param name="prevCols">The list of predecessor vertices in the POA, and their associated score vectors</param>
        /// <param name="currentVertex">The current vector in the POA that we are aligning</param>
        /// <param name="range"></param>
        /// <param name="freeInsertions">Whether to penalize for insertions or not.</param>
        /// <param name="colBestScore"> </param>
        /// <param name="colBestScoreReadPos"> </param>
        /// <param name="vertextReward">Score for matching to this node -- goes down for lower-quality nodes</param>
        /// <returns>The column of the alignment score for the current vertex</returns>
        public PoaColumn AlignCol(
            string read, Tuple<AVertex, PoaColumn>[] prevCols, AVertex currentVertex, Range range,
            bool freeInsertions, out int colBestScore, out int colBestScoreReadPos, short vertextReward)
        {
            int[] score;
            PartialOrderPath[] path;
            int offset = range.Start;

            if (buffers != null)
            {
                var bufs = buffers.GetArrays(range.Length);
                score = bufs.Item1;
                path = bufs.Item2;
            }
            else
            {
                score = new int[range.Length + 1];
                path = new PartialOrderPath[range.Length + 1];
            }

            short extra = Scores[0];
            short missing = Scores[1];
            short mismatch = Scores[2];
            //short match = Scores[3];
            short match = vertextReward;
            short branch = Scores.Length == 5 ? Scores[4] : extra;

            // If this is the last column, we let through an unbounded number of insertions, to allow partial overlaps
            if (freeInsertions)
            {
                extra = 0;
                branch = 0;
            }


            // Fill out the read==0 position.  If this vertex doesn't have any in-edges, then it's a start, otherwise it's a deletion
            if (range.Start == 0 && prevCols.Length == 0)
            {
                score[0] = 0;
                var ps = new PartialOrderPath()
                {
                    Mode = AlignMode.Start,
                    PrevVertex = null
                };
                path[0] = ps;
            }

            colBestScore = int.MinValue;
            colBestScoreReadPos = 0;

            var start = Math.Max(0, range.Start);
            var end = Math.Min(range.Stop, read.Length);

            // Align all rows in the range that we were handed
            for (int i = start + 1; i <= end; i++)
            {            
                var bestScore = 0;

                var pathSeg = new PartialOrderPath
                    {
                        Mode = AlignMode.Start,
                        PrevVertex = null
                    };

                int sc;

                var npc = prevCols.Length;
                for(int j = 0; j < npc; j++)
                {
                    var pc = prevCols[j];
                    if (pc.Item2.HasPosition(i - 1))
                    {
                        // Match 
                        var m = pc.Item1.Base == read[i - 1];
                        sc = m ? pc.Item2.Score(i - 1) + match : pc.Item2.Score(i - 1) + mismatch;

                        if (sc > bestScore)
                        {
                            bestScore = sc;
                            pathSeg.Mode = m ? AlignMode.Match : AlignMode.Mismatch;
                            pathSeg.PrevVertex = pc.Item1;
                        }
                    }

                    if (pc.Item2.HasPosition(i))
                    {

                        sc = pc.Item2.Score(i) + missing;
                        if (sc > bestScore)
                        {
                            bestScore = sc;
                            pathSeg.Mode = AlignMode.Delete;
                            pathSeg.PrevVertex = pc.Item1;
                        }
                    }
                }

                // Insert -- these come from the current segment
                var br = currentVertex.Base == read[i - 1];
                sc = score[i - 1 - offset] + (br ? branch : extra);

                if (sc > bestScore)
                {
                    bestScore = sc;
                    pathSeg.Mode = AlignMode.Insert;
                    pathSeg.PrevVertex = currentVertex;
                }

                score[i - offset] = bestScore;
                path[i - offset] = pathSeg;

                if(bestScore > colBestScore)
                {
                    colBestScore = bestScore;
                    colBestScoreReadPos = i;
                }
            }

            var col = new PoaColumn()
            {
                StartBase = range.Start,
                ScoreChunk = score,
                PathChunk = path
            };

            return col;
        }

        /*
         * Pipeline Hackers Guide
         * 
         * 
         * 1. Get Resharper, learn the keybindings
         * 2. Step through the code
         * 3. Get dotTrace and .NET memory profiler and learn to use them
         * 4. Learn to use the mono profiler
         * 5. You should know exactly what the top 20 methods in the profiles are doing
         * 6. Keep the lessons about Premature Optimization in mind. Unless your are writing code that does a significant amount of work on every pulse or frame
         *    it is not likely to show up in profiles.  Write the code to be a clear as possible. 
        */
        
        /// <summary>
        /// Find the consensus sequence of the POA, and the extents of each read with respect to that consensus
        /// </summary>
        public Dictionary<int,Tuple<AlignCell, AlignCell>> FindConsensusAndAlignments(int minCoverage, out float consensusScore, out string consensusSequence)
        {
            // A node gets a score of NumReads if all reads go through it, and a score of -NumReads if no reads go through it
            // The shift of -0.0001 breaks ties in favor of skipping half-full nodes.  In the 2 reads case this will get
            // rid of insertions which are the more common error.

            ComputeCoverage();

            Func<AVertex, float> score = v =>
            {
                if (v.DoNotCall || v.Coverage < minCoverage)
                    return Single.NegativeInfinity;

                // penalize regions with low coverage -- coverage should only be able to drop to
                // NumReads - 2 in a SMRTBell read.
                return 2 * v.NumSequences - 1 * Math.Max(NumReads-3, v.Coverage) - 0.1f;
            };

            var path = GraphAlgos.MaxPath(Graph, score);
            var seq = path.Select(v => v.Base).ToArray();

            var cScore = path.Select(v => (float)(2 * v.NumSequences - v.Coverage) / v.Coverage).Average();
            consensusScore = cScore;

            consensusSequence = new string(seq);

            var startPos = new Dictionary<int, Tuple<int, int>>();
            var endPos = new Dictionary<int, Tuple<int, int>>();

            for(int i = 0; i < path.Count; i++)
            {
                var node = path[i];

                foreach (var rp in node.ReadPointers)
                {
                    if (!startPos.ContainsKey(rp.ReadId))
                    {
                        startPos[rp.ReadId] = new Tuple<int, int>(rp.ReadPosition, i);
                    }
                }
            }

            for (int i = path.Count - 1; i >= 0; i--)
            {
                var node = path[i];

                foreach (var rp in node.ReadPointers)
                {
                    if (!endPos.ContainsKey(rp.ReadId))
                    {
                        endPos[rp.ReadId] = new Tuple<int, int>(rp.ReadPosition, i);
                    }
                }
            }

            var extentDictionary = new Dictionary<int, Tuple<AlignCell, AlignCell>>();

            foreach(var i in Enumerable.Range(0, NumReads))
            {
                if (startPos.ContainsKey(i))
                {
                    var start = startPos[i];
                    var end = endPos[i];

                    var extents = new Tuple<AlignCell, AlignCell>(
                        new AlignCell(AlignMode.Match, start.Item1, start.Item2),
                        new AlignCell(AlignMode.Match, end.Item1, end.Item2));

                    extentDictionary[i] = extents;
                }
            }

            return extentDictionary;
        }


        /// <summary>
        /// Write the POA to a dot file, suitable for visualization with Graphviz
        /// </summary>
        /// <param name="filename"></param>
        public void WriteGraphViz(string filename)
        {
            var graphviz = new GraphvizAlgorithm<AVertex, Edge<AVertex>>(Graph);
            graphviz.FormatVertex += Format;

            // render
            string output = graphviz.Generate();
            System.IO.File.WriteAllText(filename, output);
        }

        
        public void Format(Object sender, FormatVertexEventArgs<AVertex> e)
        {
            e.VertexFormatter.Label = e.Vertex.ToString();
        }

        public struct BitSet
        {
            private long bits;

            private BitSet(long bits)
            {
                this.bits = bits;
            }

            public BitSet(IEnumerable<ushort> v)
            {
                bits = 0;

                foreach (var i in v)
                {
                    bits = bits | (1L << i);
                }
            }
            
            public BitSet Add(ushort item)
            {
                return new BitSet(bits | (1L << item));
            }


            public BitSet UnionWith(BitSet other)
            {
                return new BitSet(bits | other.bits);
            }

            public BitSet IntersectWith(BitSet other)
            {
                return new BitSet(bits & other.bits);
            }


            public IEnumerable<ushort> Items
            {
                get
                {
                    for (int i = 0; i < 64; i++)
                    {
                        if (((bits >> i) & 1L) == 1L)
                        {
                            yield return (ushort) i;
                        }
                    }
                }
            }

            public bool Contains(ushort i)
            {
                if (((bits >> i) & 1L) == 1L)
                {
                    return true;
                }
                return false;
            }

            public int Count
            {
                get
                {
                    int count = 0;
                    for (int i = 0; i < 64; i++)
                    {
                        if (((bits >> i) & 1L) == 1L)
                        {
                            count++;
                        }
                    }
                    return count;
                }
            }

        }


        /// <summary>
        /// Annotate the POA with local coverage information.  If a read 'passes-by' a node, then it counts toward the 
        /// coverage of that node.
        /// </summary>
        public void ComputeCoverage()
        {
            // Array of AVertex objects
            var tNodes = Graph.TopologicalSort;

            // Lookup from a vertex to it's position in tNods
            var vertexDict = Graph.VertexDict;

            var incomingSets = new BitSet[tNodes.Length];
            var outgoingSets = new BitSet[tNodes.Length];

            Func<AVertex, BitSet> incoming = vtx => incomingSets[vertexDict[vtx]];
            Func<AVertex, BitSet> outgoing = vtx => outgoingSets[vertexDict[vtx]];

            var readsSeen = new BitSet();

            // Figure out the set of reads that start at a predecessor of this node  
            foreach (var v in tNodes)
            {
                var newReads = v.ReadPointers.Where(rp => !readsSeen.Contains(rp.ReadId)).Select(rp => rp.ReadId);
                var incomingSet = new BitSet(newReads);
                readsSeen = readsSeen.UnionWith(incomingSet);

                foreach(var inBound in Graph.InEdges(v).Select(e => incoming(e.Source)))
                {
                    incomingSet = incomingSet.UnionWith(inBound);
                }

                incomingSets[vertexDict[v]] = incomingSet;
            }

            readsSeen = new BitSet();

            // Do the same procedure as above, except in the reverse direction
            for (int i = tNodes.Length - 1; i >= 0; i--)
            {
                var v = tNodes[i];
                var newReads = v.ReadPointers.Where(rp => !readsSeen.Contains(rp.ReadId)).Select(rp => rp.ReadId);
                var outgoingSet = new BitSet(newReads);
                readsSeen = readsSeen.UnionWith(outgoingSet);

                foreach (var outBound in Graph.OutEdges(v).Select(e => outgoing(e.Target)))
                {
                    outgoingSet = outgoingSet.UnionWith(outBound);
                }
                outgoingSets[vertexDict[v]] = outgoingSet;
            }

            // Compute the 'coverage' for each node.  The coverage is the number of reads in the
            // intersection of the incoming and reverse incoming sets.
            for(int i = 0; i < tNodes.Length; i++)
            {
                var common = incomingSets[i].IntersectWith(outgoingSets[i]);
                tNodes[i].Coverage = common.Count;
            }
        }

        public void CleanPoa(int minCoverage)
        {
            Graph.RemoveVertexIf(v => v.Coverage == 1 || (v.NumSequences == 1 & v.Coverage > minCoverage));
        }



        /// <summary>
        /// Structure for tracking the banding extents of POA alignment
        /// </summary>
        public struct Range
        {
            public Range(int start, int stop)
            {
                Start = start;
                Stop = stop;
            }

            public bool Equals(Range other)
            {
                return Start == other.Start && Stop == other.Stop;
            }

            public override bool Equals(object obj)
            {
                if (ReferenceEquals(null, obj))
                {
                    return false;
                }
                return obj is Range && Equals((Range)obj);
            }

            public override int GetHashCode()
            {
                unchecked
                {
                    return (Start * 397) ^ Stop;
                }
            }

            public readonly int Start;
            public readonly int Stop;

            public int Length
            {
                get { return Math.Max(0, Stop - Start); }
            }

            public Range UnionWith(Range r2)
            {
                var r1 = this;
                return new Range(Math.Min(r1.Start, r2.Start), Math.Max(r1.Stop, r2.Stop));
            }

            public Range Next(int readLength, int squeeze = 0)
            {
                if (this == Empty)
                    return Empty;

                return new Range(Math.Min(Start + 1 + squeeze, readLength), Math.Min(Stop + 1 - squeeze, readLength));
            }

            public Range Prev(int readLength, int squeeze = 0)
            {
                if (this == Empty)
                    return Empty;

                return new Range(Math.Max(0, Start - 1 + squeeze), Math.Max(0, Stop - 1 - squeeze));
            }

            public static Range Empty
            {
                get { return new Range(Int32.MaxValue/2, Int32.MinValue/2 ); }
            }

            public static Range Union(IEnumerable<Range> ranges)
            {
                return ranges.Aggregate(Empty, (r1, r2) => r1.UnionWith(r2));
            }

            public static Range Union(Range r1, Range r2)
            {
                return new Range(Math.Min(r1.Start, r2.Start), Math.Max(r1.Stop, r2.Stop));
            }

            public static bool operator ==(Range r1, Range r2)
            {
                return r1.Start == r2.Start && r1.Stop == r2.Stop;
            }

            public static bool operator !=(Range r1, Range r2)
            {
                return !(r1 == r2);
            }
        }

        /// <summary>
        /// Perform sparse alignments from fragments of new reads against reads already added to the POA, and
        /// use them to set the band to use during full POA alignment.
        /// </summary>
        public class AlignmentBackbone
        {
            public static int Width = 30;
            private string[] reads;

            /// <summary>
            /// Stores the list of reads with which we will form the POA
            /// </summary>
            /// <param name="reads"></param>
            public AlignmentBackbone(string[] reads)
            {
                this.reads = reads;
            }

            /// <summary>
            /// Find a fragment in fragments whose Y position matches readPosition, and return it,
            /// otherwise return null
            /// </summary>
            public static Fragment FindFragment(int readPosition, Fragment[] fragments)
            {
                // Do a binary search for a fragment with the given readPosition in the Y field
                int start = 0;
                int end = fragments.Length - 1;

                while (end - start > 1)
                {
                    var t = (end + start) / 2;
                    if (fragments[t].Y == readPosition)
                    {
                        return fragments[t];
                    }
                    else if (fragments[t].Y < readPosition)
                    {
                        start = t;
                    }
                    else
                    {
                        end = t;
                    }
                }

                return null;
            }

            /// <summary>
            /// Call before the read alignment to get a closure that supplies convert POA nodes
            /// (indexed by topo sort) into the range of the read that should be aligned.
            /// </summary>
            public Func<int, Range> PrepareRanges(int readId,  PoaGraph<AVertex, Edge<AVertex>> graph)
            {
                var topo = graph.TopologicalSort;
                var vertexDict = graph.VertexDict;
                var readLength = reads[readId].Length;

                // Get fragment alignments of this read with all prior reads
                var sparseAl = new SparseAligner(6);
                Fragment[][] frags = (readId).Fill(idx => sparseAl.SparseAlignStrings(reads[readId], reads[idx]).ToArray());

                // Find any fragments that align the current read to bases in the ith node of the POA topo sort.
                // Return a range that surrounds these fragments.
                Func<int, Range?> findDirectRange = i =>
                {
                    // see if any of the fragment alignments have fragments corresponding to the bases aligned to this POA node
                    var readPointers = topo[i].ReadPointers;

                    var directRange =
                        readPointers.Select(rp =>
                        {
                            var fragArray = frags[rp.ReadId];
                            var matchFrag = FindFragment(rp.ReadPosition, fragArray);

                            if (matchFrag != null)
                            {
                                return new Range(Math.Max(0, matchFrag.X - Width), Math.Min(reads[readId].Length, matchFrag.X + Width));
                            }
                            else
                            {
                                return Range.Empty;
                            }
                        }).Aggregate(Range.Empty, Range.Union);

                    if (directRange == Range.Empty)
                    {
                        return null;
                    }
                    else
                    {
                        return directRange;
                    }

                };

                var fwdMarks = new Range[topo.Length];
                var revMarks = new Range[topo.Length];

                // We perform a forward and backward sweep of the POA to set the bounds of the POA alignment
                // The forward (backward) ranges are set as follows:
                // 1. If the current node has fragment hits, set the bounds to the union of these ranges
                // 2. If the current node doesn't fragment hits, it's the union of the stepped forward (backward) predecessor (successor) hits
                // The range for a node is the union of the forward and backward ranges set in this manner.

                for (int i = 0; i < fwdMarks.Length; i++)
                {
                    var freshRange = findDirectRange(i);

                    if (freshRange.HasValue)
                    {
                        // If we have a fresh range, use it
                        fwdMarks[i] = freshRange.Value;
                    }
                    else
                    {
                        // Otherwise propagate the previous ranges
                        var inNodes = graph.InEdges(topo[i]).Select(e => e.Source);
                        var inRanges = inNodes.Select(v => fwdMarks[vertexDict[v]].Next(readLength));
                        var range = Range.Union(inRanges);
                        fwdMarks[i] = range;
                    }
                }


                for (int i = fwdMarks.Length - 1; i >= 0; i--)
                {
                    var freshRange = findDirectRange(i);

                    if (freshRange.HasValue)
                    {
                        // If we have a fresh range, use it
                        revMarks[i] = freshRange.Value;
                    }
                    else
                    {
                        // Otherwise propagate the previous ranges
                        var outNodes = graph.OutEdges(topo[i]).Select(e => e.Target);
                        var outRanges = outNodes.Select(v => revMarks[vertexDict[v]].Prev(readLength));
                        var range = Range.Union(outRanges);
                        revMarks[i] = range;
                    }
                }

                return i => fwdMarks[i].UnionWith(revMarks[i]);
            }

            /// <summary>
            /// Call before the read alignment to get a closure that supplies convert POA nodes
            /// (indexed by topo sort) into the range of the read that should be aligned.
            /// </summary>
            public static Func<int, Range> SimpleRange(string read, string consensus, PoaGraph<AVertex, Edge<AVertex>> graph, AVertex[] consensusPath)
            {
                var topo = graph.TopologicalSort;
                var vertexDict = graph.VertexDict;

                var readLength = read.Length;

                // Get fragment alignments of this read with all prior reads
                var sparseAl = new SparseAligner(6);
                Fragment[] frags = sparseAl.SparseAlignStrings(read, consensus).ToArray();
                
                var consensusPosition = new Dictionary<AVertex, int>();
                consensusPath.ForEach((idx, vtx) => consensusPosition[vtx] = idx);


                // Find any fragments that align the current read to bases in the ith node of the POA topo sort.
                // Return a range that surrounds these fragments.
                Func<int, Range?> findDirectRange = i =>
                    {
                        int pos;
                        
                        if(consensusPosition.TryGetValue(topo[i], out pos))
                        {
                            var matchFrag = FindFragment(pos, frags);
                            
                            if (matchFrag != null)
                            {
                                return new Range(Math.Max(0, matchFrag.X - Width), Math.Min(read.Length, matchFrag.X + Width));
                            }
                            else
                            {
                                return null;
                            }
                        }

                        return null;
                    };


                var fwdMarks = new Range[topo.Length];
                var revMarks = new Range[topo.Length];

                // We perform a forward and backward sweep of the POA to set the bounds of the POA alignment
                // The forward (backward) ranges are set as follows:
                // 1. If the current node has fragment hits, set the bounds to the union of these ranges
                // 2. If the current node doesn't fragment hits, it's the union of the stepped forward (backward) predecessor (successor) hits
                // The range for a node is the union of the forward and backward ranges set in this manner.

                for (int i = 0; i < fwdMarks.Length; i++)
                {
                    var freshRange = findDirectRange(i);

                    if (freshRange.HasValue)
                    {
                        // If we have a fresh range, use it
                        fwdMarks[i] = freshRange.Value;
                    }
                    else
                    {
                        // Otherwise propagate the previous ranges
                        var inRanges = graph.InEdges(topo[i]).Select(e => fwdMarks[vertexDict[e.Source]].Next(readLength, 1));
                        var range = Range.Union(inRanges);
                        fwdMarks[i] = range;
                    }
                }


                for (int i = fwdMarks.Length - 1; i >= 0; i--)  
                {
                    var freshRange = findDirectRange(i);

                    if (freshRange.HasValue)
                    {
                        // If we have a fresh range, use it
                        revMarks[i] = freshRange.Value;
                    }
                    else
                    {
                        // Otherwise propagate the previous ranges
                        var outRanges = graph.OutEdges(topo[i]).Select(e => revMarks[vertexDict[e.Target]].Prev(readLength, 1));
                        var range = Range.Union(outRanges);
                        revMarks[i] = range;
                    }
                }

                return i => fwdMarks[i].UnionWith(revMarks[i]);
            }

        }
    }

    public struct ReadPointer
    {
        public bool ReadStart { get; set; }
        public bool ReadEnd { get; set; }
        public ushort ReadId;
        public ushort ReadPosition;
    }

    public class AVertex
    {
        // DNA base of this vertex
        public char Base;

        // Local coverage of reads passing this vertex
        public int Coverage;

        // Do not include this base in a consensus sequence
        public bool DoNotCall = false;
        
        // Indicates which bases from which reads are mapped into this vertex
        public List<ReadPointer> ReadPointers = new List<ReadPointer>();

        // Number of reads mapping into this node ( == ReadPointers.Length)
        public ushort NumSequences = 0;

        public AVertex(char bse)
        {
            Base = bse;
            NumSequences = 0;
        }

        public AVertex(char bse, int readId, int readPos) : this(bse)
        {
            AddSequence(readId, readPos);
        }

        public void AddSequence(int readId, int readPos)
        {
            ReadPointers.Add(new ReadPointer { ReadId =  (ushort) readId, ReadPosition = (ushort) readPos });
            NumSequences++;
        }

        public new string ToString()
        {
            if (DoNotCall)
                return String.Format("{0} n={1} c={2} -- NO CALL", Base, NumSequences, Coverage);

            return String.Format("{0} n={1} c={2}", Base, NumSequences, Coverage);
        }
    }

    public class PoaGraph<TNode, TEdge> : BidirectionalGraph<TNode, TEdge> where TEdge : IEdge<TNode>
    {
        private TNode[] topologicalSort;
        private Dictionary<TNode, int> vertexDict; 

        public PoaGraph(bool allowParallelEdges) : base(allowParallelEdges)
        {
            
        }

        private void ClearTopo()
        {
            topologicalSort = null;
            vertexDict = null;
        }

        private void ComputeTopo()
        {
            var topo = new C5.LinkedList<TNode>();
            this.TopologicalSort(topo);

            topologicalSort = topo.ToArray();

            vertexDict = new Dictionary<TNode, int>();
            topologicalSort.ForEach((idx, v) => vertexDict.Add(v, idx));
        }

        public TNode[] TopologicalSort
        {
            get
            {
                if (topologicalSort == null)
                    ComputeTopo();

                return topologicalSort;
            }
        }

        public Dictionary<TNode, int> VertexDict
        {
            get
            {
                if(vertexDict == null)
                    ComputeTopo();

                return vertexDict;
            }
        }

        protected override void OnVertexAdded(TNode args)
        {
            base.OnVertexAdded(args);
            ClearTopo();
        }

        protected override void OnVertexRemoved(TNode args)
        {
            base.OnVertexRemoved(args);
            ClearTopo();
        }


        protected override void OnEdgeRemoved(TEdge args)
        {
            base.OnEdgeRemoved(args);
            ClearTopo();
        }
    }



    public class ChunkVertex<TVertex> where TVertex : class
    {
        public List<TVertex> Chunk;
        public ChunkVertex(List<TVertex> vertices)
        {
            Chunk = vertices;
            if(vertices.Count == 0 && System.Diagnostics.Debugger.IsAttached)
                System.Diagnostics.Debugger.Break();

        }

        public bool TrySplit(BidirectionalGraph<ChunkVertex<TVertex>, IEdge<ChunkVertex<TVertex>>> graph, TVertex splitVertex, List<TVertex> newChunk, bool outBound)
        {
            // Scan my chunk to see if we have the vertx
            int index = Chunk.FindIndex(v => v == splitVertex);
            if(index >= 0)
            {
                // Tailing edge
                if (outBound)
                {
                    var outEdges = graph.OutEdges(this).ToArray();

                    var myNewChunk = Chunk.Take(index + 1).ToList();
                    var descChunk = Chunk.Skip(index + 1).ToList();

                    // Give myself the shortened chunk
                    Chunk = myNewChunk;

                    // Make the child chunk, and transfer all my children to it
                    var splitChildVertex = new ChunkVertex<TVertex>(descChunk);
                    graph.AddVertex(splitChildVertex);

                    foreach (var outEdge in outEdges)
                    {
                        var e = new Edge<ChunkVertex<TVertex>>(splitChildVertex, outEdge.Target);
                        graph.AddEdge(e);
                        graph.RemoveEdge(outEdge);
                    }


                    // Make the new child vertex
                    var newChildVertex = new ChunkVertex<TVertex>(newChunk);
                    graph.AddVertex(newChildVertex);

                    graph.AddEdge(new Edge<ChunkVertex<TVertex>>(this, splitChildVertex));
                    graph.AddEdge(new Edge<ChunkVertex<TVertex>>(this, newChildVertex));
                }
                else
                {
                    var inEdges = graph.InEdges(this).ToArray();

                    var predChunk = Chunk.Take(index).ToList();
                    var myNewChunk = Chunk.Skip(index).ToList();

                    // Give myself the shortened chunk
                    Chunk = myNewChunk;

                    // Make the child chunk, and transfer all my children to it
                    var splitParentVertex = new ChunkVertex<TVertex>(predChunk);
                    graph.AddVertex(splitParentVertex);

                    foreach (var inEdge in inEdges)
                    {
                        var e = new Edge<ChunkVertex<TVertex>>(inEdge.Source, splitParentVertex);
                        graph.AddEdge(e);
                        graph.RemoveEdge(inEdge);
                    }


                    // Make the new child vertex
                    var newParentVertex = new ChunkVertex<TVertex>(newChunk);
                    graph.AddVertex(newParentVertex);

                    graph.AddEdge(new Edge<ChunkVertex<TVertex>>(splitParentVertex, this));
                    graph.AddEdge(new Edge<ChunkVertex<TVertex>>(newParentVertex, this));
                }

                return true;
            }
            else
            {
                return false;
            }
        }

        public new string ToString()
        {
            var t = typeof (TVertex);

            if(t == typeof(AVertex))
            {
                var chunk = Chunk as List<AVertex>;

                var seq = chunk.Select(v => v.Base).ToArray();
                var score = chunk.Select(v => 2.0*v.NumSequences - v.Coverage).Average();

                return String.Format("{0} \r\n Score: {1}", seq, score);
            }

            return "";

        }

    }
}
