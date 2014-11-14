using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using PacBio.Align;
using PacBio.IO;
using PacBio.Utils;
using QuickGraph;
using QuickGraph.Algorithms;
//using QuickGraph.Graphviz;

namespace PacBio.Consensus
{
    using MutationType = ConsensusCore.MutationType;
    using PartialOrderPath = PacBio.Consensus.PartialOrderPath<PacBio.Consensus.Vertex>;

    public static class GraphAlgos
    {
        /// <summary>
        /// Find the highest scoring path that exists in a DAG. The scoring function provides a (possibly negative) score for including a
        /// vertex in the path.
        /// </summary>
        /// <param name="graph">A DAG</param>
        /// <param name="f">A function giving the score of a vertex</param>
        /// <returns>The list of vertexs in the maximum scoring path, in order</returns>
        public static List<TVertex> MaxPath<TVertex, TEdge>(BidirectionalGraph<TVertex, TEdge> graph, Func<TVertex, float> f) where TEdge : IEdge<TVertex>
        {
            var scores = new float[graph.VertexCount];

            var path = Enumerable.Repeat(-1, graph.VertexCount).ToArray();

            TVertex[] vertices;
            Dictionary<TVertex, int> dictIdx;

            if (graph is PoaGraph<TVertex, TEdge>)
            {
                var g = graph as PoaGraph<TVertex, TEdge>;
                vertices = g.TopologicalSort;
                dictIdx = g.VertexDict;
            }
            else
            {
                vertices = graph.TopologicalSort().ToArray();
                dictIdx = new Dictionary<TVertex, int>();
            }

            for(int i =0; i < vertices.Length; i++)
                dictIdx[vertices[i]] = i;

            for(int i = 0; i < vertices.Length; i++)
            {
                var vertex = vertices[i];
                var vertexScore = f(vertex);

                var inEdges = graph.InEdges(vertex).ToArray();

                scores[i] = vertexScore;
                path[i] = -1;

                for(int j = 0; j < inEdges.Length; j++)
                {
                    var src = inEdges[j].Source;
                    var srcIdx = dictIdx[src];

                    var sc = scores[srcIdx] + vertexScore;
                    if (sc > scores[i])
                    {
                        scores[i] = sc;
                        path[i] = srcIdx;
                    }
                }
            }

            // Run the trace back to get the best score
            var pathVertices = new C5.LinkedList<TVertex>();
            var nextIdx = scores.IMax();

            while(nextIdx >= 0)
            {
                pathVertices.Insert(0, vertices[nextIdx]);
                nextIdx = path[nextIdx];
            }

            return pathVertices.ToList();
        }



        /// <summary>
        /// Find the highest scoring path that exists in a DAG. The scoring function provides a (possibly negative) score for including a
        /// vertex in the path.
        /// </summary>
        /// <param name="graph">A DAG</param>
        /// <param name="f">A function giving the score of a vertex</param>
        /// <returns>The list of vertexs in the maximum scoring path, in order</returns>
        public static BidirectionalGraph<ChunkVertex<TVertex>, IEdge<ChunkVertex<TVertex>>> MaxMultiPath<TVertex, TEdge>(BidirectionalGraph<TVertex, TEdge> graph, Func<TVertex, float> f)
            where TEdge : IEdge<TVertex> 
            where TVertex : class
        {
            var scores = new float[graph.VertexCount];
            var path = Enumerable.Repeat(-1, graph.VertexCount).ToArray();

            var vertexInPath = new bool[graph.VertexCount];
            var predInPath = new bool[graph.VertexCount];
            var descInPath = new bool[graph.VertexCount];

            var vertices = graph.TopologicalSort().ToArray();
            var dictIdx = new Dictionary<TVertex, int>();

            for (int i = 0; i < vertices.Length; i++)
                dictIdx[vertices[i]] = i;

            for (int i = 0; i < vertices.Length; i++)
            {
                var vertex = vertices[i];
                var vertexScore = f(vertex);

                var inEdges = graph.InEdges(vertex).ToArray();

                scores[i] = vertexScore;
                path[i] = -1;

                for (int j = 0; j < inEdges.Length; j++)
                {
                    var src = inEdges[j].Source;
                    var srcIdx = dictIdx[src];

                    var sc = scores[srcIdx] + vertexScore;
                    if (sc > scores[i])
                    {
                        scores[i] = sc;
                        path[i] = srcIdx;
                    }
                }
            }

            var chunkGraph = new BidirectionalGraph<ChunkVertex<TVertex>, IEdge<ChunkVertex<TVertex>>>();

            // Decode all the paths.
            while (true)
            {
                // Find the best path
                var pathVertices = new List<TVertex>();

                var nextIdx = -1;
                var bestScore = 0.0f;

                // Find the best score that hasn't already been found in a path.
                for (int i = 0; i < vertices.Length; i++)
                {
                    if((!predInPath[i] || !descInPath[i]) && scores[i] > bestScore)
                    {
                        nextIdx = i;
                        bestScore = scores[i];
                    }
                }

                // Can't find any new paths to take
                if (nextIdx == -1)
                    break;

                // Trace back until we get to the end of the path, or we merge into another path
                while (nextIdx >= 0 && !vertexInPath[nextIdx])
                {
                    vertexInPath[nextIdx] = true;
                    predInPath[nextIdx] = true;
                    descInPath[nextIdx] = true;
                    pathVertices.Insert(0, vertices[nextIdx]);
                    nextIdx = path[nextIdx];
                }

                if (nextIdx == -1)
                {
                    if (chunkGraph.VertexCount == 0)
                    {
                        var mainChunk = new ChunkVertex<TVertex>(pathVertices);
                        chunkGraph.AddVertex(mainChunk);
                    }
                    // Don't handle front-side spitting right now
                }
                else if(bestScore - scores[nextIdx] > 0)
                {
                    var foundSplit = chunkGraph.Vertices.Any(vertex => vertex.TrySplit(chunkGraph, vertices[nextIdx], pathVertices, true));

                    if (foundSplit == false)
                        throw new Exception("couldn't split chunk correctly");
                }


                // update the path-taken markers -- we can't start a new path from a node that
                // has been 'covered' by a previous path
                for(int i = 0; i < vertices.Length; i++)
                {
                    var vertex = vertices[i];
                    var inPath = predInPath[i];

                    if (!inPath)
                    {
                        var inEdges = graph.InEdges(vertex);
                        if (inEdges.Any(inEdge => predInPath[dictIdx[inEdge.Source]]))
                        {
                            descInPath[i] = true;
                        }
                    }
                }

                for (int i = vertices.Length - 1; i >= 0; i--)
                {
                    var vertex = vertices[i];
                    var inPath = descInPath[i];

                    if (!inPath)
                    {
                        var outEdges = graph.OutEdges(vertex);
                        if (outEdges.Any(outEdge => descInPath[dictIdx[outEdge.Target]]))
                        {
                            predInPath[i] = true;
                        }
                    }
                }
            }

            return chunkGraph;
        }
     }

    public struct PartialOrderPath<T>
    {
        public AlignMode Mode;
        public T PrevVertex;

        public static PartialOrderPath<T> Default
        {
            get
            {
                return new PartialOrderPath<T>()
                {
                    Mode = AlignMode.Start,
                    PrevVertex = default(T)
                };
            }
        }

    }
    
    public class Vertex
    {
        public char Base;
        public int NumSequences;
        public int Coverage;

        public List<int> startReads;
        public List<int> endReads;

        public Vertex(char bse)
        {
            Base = bse;
            NumSequences = 1;
        }

        public void AddSequence(int s)
        {
            NumSequences++;
        }

        public void ReadStart(int read)
        {
            if(startReads == null)
            {
                startReads = new List<int>();
            }

            startReads.Add(read);
        }

        public void ReadEnd(int read)
        {
            if (endReads == null)
            {
                endReads = new List<int>();
            }

            endReads.Add(read);
        }

        public new string ToString()
        {
            return String.Format("{0} - n={1}, cov={2}", Base, NumSequences, Coverage);
        }
    }
}
