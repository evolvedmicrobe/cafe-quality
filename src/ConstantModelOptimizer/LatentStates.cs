using System;

namespace ConstantModelOptimizer
{
    /// <summary>
    /// All possible latent states, so we can keep them in one matrix
    /// </summary>
    public struct LatentStates
    {
        /// <summary>
        /// Always in log space.
        /// </summary>
        public double Match, Stick, Branch, Dark, Merge, Total;
        public LatentStates()
        {
            Match = Double.NegativeInfinity;
            Stick = Double.NegativeInfinity;
            Branch = Double.NegativeInfinity;
            Dark = Double.NegativeInfinity;
            Merge = Double.NegativeInfinity;
            Total = Double.NegativeInfinity;
        }
        public void Clear() {
            Match = Double.NegativeInfinity;
            Stick = Double.NegativeInfinity;
            Branch = Double.NegativeInfinity;
            Dark = Double.NegativeInfinity;
            Merge = Double.NegativeInfinity;
            Total = Double.NegativeInfinity;
        }
        public void SetTotal()
        {
            Total = MathUtils.logsumlog (Match, Stick, Branch, Dark, Merge);
        }

        /// <summary>
        /// Combines element wise to produce the log of the exp sum, e.g.
        /// Match = log( exp(Match_{1}) + exp(Match_{2}))
        /// The total is cleared
        /// </summary>
        /// <param name="other">Other.</param>
        public void Addin(LatentStates other)
        {
            Match = MathUtils.logsumlog (Match, other.Match);
            Merge = MathUtils.logsumlog (Merge, other.Merge);
            Branch = MathUtils.logsumlog (Branch, other.Branch);
            Stick = MathUtils.logsumlog (Stick, other.Stick);
            Dark = MathUtils.logsumlog (Dark, other.Dark);
            Total = Double.NegativeInfinity;
        }

        /// <summary>
        /// This subtracts a constant from all states, 
        /// in practice used for normalizing by the total probability
        /// in Bayes rule.
        /// </summary>
        /// <param name="val">Value.</param>
        public void RemoveConstant(double val) {
            Match -= val;
            Merge -= val;
            Branch -= val;
            Stick -= val;
            Dark -= val;
            Total = Double.NegativeInfinity;
        }
    }
}

