using System;
using System.Collections.Generic;
using System.Linq;

namespace Bio.Algorithms.Assembly.Graph
{
    /// <summary>
    /// Represents a path in De Bruijn graph.
    /// </summary>
    public class DeBruijnPath
    {
        /// <summary>
        /// List of node in De Bruijn graph path.
        /// </summary>
        private List<DeBruijnNode> path;

        /// <summary>
        /// Initializes a new instance of the DeBruijnPath class.
        /// </summary>
        public DeBruijnPath()
        {
            this.path = new List<DeBruijnNode>();
        }

        /// <summary>
        /// Initializes a new instance of the DeBruijnPath class with specified nodes.
        /// </summary>
        /// <param name="nodes">List of nodes.</param>
        public DeBruijnPath(IEnumerable<DeBruijnNode> nodes)
        {
            this.path = new List<DeBruijnNode>(nodes);
        }

        /// <summary>
        /// Initializes a new instance of the DeBruijnPath class with specified node.
        /// </summary>
        /// <param name="node">Graph node.</param>
        public DeBruijnPath(DeBruijnNode node)
        {
            this.path = new List<DeBruijnNode> { node };
        }

        /// <summary>
        /// Gets the list of nodes in path.
        /// </summary>
        public IList<DeBruijnNode> PathNodes
        {
            get { return this.path; }
        }

        /// <summary>
        /// Removes all nodes from path that match the given predicate.
        /// </summary>
        /// <param name="predicate">Predicate to remove nodes.</param>
        public void RemoveAll(Predicate<DeBruijnNode> predicate)
        {
            this.path.RemoveAll(predicate);
        }

        /// <summary>
        /// Returns the original sequence by moving through the path.
        /// 
        /// Assumes all paths are connected and that we know which way to go.
        /// </summary>
        /// <returns></returns>
        public Sequence ConvertToSequence(int kmerLength)
        {            
            List<byte> sequence = new List<byte>(this.PathNodes.Count + kmerLength - 1);
            if (PathNodes.Count == 0)
            {
                return null; ;
            }
            else if (PathNodes.Count == 1)
            {
                return new Sequence(Alphabets.DNA, PathNodes[0].GetOriginalSymbols(kmerLength));
            }
            else
            {
                //  Build up the sequence by adding nodes.

                //  For the first node we add its sequence, and the neighbor could be on the left or right side
                var cur_node = PathNodes[0];
                var leftNodes = cur_node.GetLeftExtensionNodesWithOrientation().Where(x => x.Key == PathNodes[1]).ToList();
                var rightNodes = cur_node.GetRightExtensionNodesWithOrientation().Where(x => x.Key == PathNodes[1]).ToList();
                var goingLeft = leftNodes.Count == 1;
                if (!((leftNodes.Count == 1) ^ (rightNodes.Count == 1)))
                {
                    throw new InvalidProgramException();
                }
                if (goingLeft)
                {
                    sequence.AddRange(cur_node.GetReverseComplementOfOriginalSymbols(kmerLength));
                }
                else
                {
                    sequence.AddRange(cur_node.GetOriginalSymbols(kmerLength));
                }
                var sameOrientation = goingLeft ? leftNodes.First().Value : rightNodes.First().Value;
                /* After the first node, we know if the correct node is on the right or left side.
                   going to loop through and add sequences */
                for(int i = 1; i < PathNodes.Count; i++)
                {
                    var grab_last_base = goingLeft ^ sameOrientation;
                    var next_node = PathNodes[i];
                   
                    // if we aren't done yet, decide which way we go next.
                    if (i < (PathNodes.Count - 1))
                    {
                        byte nextSymbol = GetNextSymbol(next_node, kmerLength, !grab_last_base);
                        sequence.Add(nextSymbol);
                        goingLeft = !grab_last_base;
                        var nextNodes = goingLeft ? next_node.GetLeftExtensionNodesWithOrientation() :
                                                     next_node.GetRightExtensionNodesWithOrientation();
                        var nextNode = nextNodes.Where(x => x.Key == PathNodes[i + 1]).FirstOrDefault();
                        if (nextNode.Key == null)
                        {
                            throw new InvalidProgramException();
                        }
                        sameOrientation = nextNode.Value;
                    }
                    else
                    {
                        //add everything for the last base
                        var bytes = grab_last_base ? next_node.GetOriginalSymbols(kmerLength) : next_node.GetReverseComplementOfOriginalSymbols(kmerLength);
                    }
                }
                return new Sequence(Alphabets.DNA, sequence.ToArray());                
            }
        }

	    /// <summary>
        /// This gets the next symbol from a node while forming chains.  This can be made a lot more efficient if it turns in to a bottleneck.
        /// all chains are extended from either the first or last base present in the node, and this base is either forward
        /// or reverse complimented, this method reflects this.
        /// </summary>
        /// <param name="node">Next node</param>
        /// <param name="graph">Graph to get symbol from</param>
        /// <param name="GetFirstNotLast">First or last base?</param>
        /// <param name="ReverseComplimentBase">Should the compliment of the base be returned</param>
        /// <returns></returns>
        private static byte GetNextSymbol(DeBruijnNode node, int kmerLength, bool GetRCofFirstBaseInsteadOfLastBase)
        {
            if (node == null)
            { throw new ArgumentNullException("node"); }
            byte[] symbols = node.GetOriginalSymbols(kmerLength);
            byte value = GetRCofFirstBaseInsteadOfLastBase ? symbols.First() : symbols.Last();
            if (GetRCofFirstBaseInsteadOfLastBase)
            {
                byte value2;
                bool rced = DnaAlphabet.Instance.TryGetComplementSymbol(value, out value2);
                //Should never happend
                if (!rced)
                {
                    throw new Exception("Could not revcomp base during graph construction");
                }
                value = value2;
            }
            return value;
        }

    }
}
