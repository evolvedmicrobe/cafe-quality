using System;

namespace ConstantModelOptimizer
{
    /// <summary>
    /// All possible latent states, so we can keep them in one matrix
    /// </summary>
    public struct LatentStates
    {
        public double Match, Stick, Branch, Dark, Merge, Total;
        public LatentStates()
        {
            Clear ();

        }
        public void Clear() {
            Match = Double.NegativeInfinity;
            Stick = Double.NegativeInfinity;
            Branch = Double.NegativeInfinity;
            Dark = Double.NegativeInfinity;
            Merge = Double.NegativeInfinity;
            Total = Double.NegativeInfinity;
        }
    }
}

