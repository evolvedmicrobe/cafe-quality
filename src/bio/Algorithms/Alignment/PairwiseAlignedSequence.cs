using System.Linq;
using System.Text;

namespace Bio.Algorithms.Alignment
{
    /// <summary>
    /// PairwiseAlignedSequence is a class containing the single aligned unit of pairwise alignment.
    /// </summary>
    public class PairwiseAlignedSequence : AlignedSequence
    {
        #region Private constants
        /// <summary>
        /// Constant string indicating consensus in meta-data.
        /// </summary>
        private const string ConsensusKey = "Consensus";

        /// <summary>
        /// Constant string indicating alignment score in meta-data.
        /// </summary>
        private const string ScoreKey = "Score";

        /// <summary>
        /// Constant string indicating offset of first sequence in alignment.
        /// </summary>
        private const string FirstOffsetKey = "FirstOffset";

        /// <summary>
        /// Constant string indicating offset of second sequence in alignment.
        /// </summary>
        private const string SecondOffsetKey = "SecondOffset";
        #endregion

        #region Constructors

        /// <summary>
        /// Initializes a new instance of the PairwiseAlignedSequence class.
        /// </summary>
        public PairwiseAlignedSequence()
            : base()
        {
            // No impl.
        }

        /// <summary>
        /// Initializes a new instance of the PairwiseAlignedSequence class
        /// Internal constructor for creating new instance of 
        /// PairwiseAlignedSequence from specified IAlignedSequence.
        /// </summary>
        /// <param name="alignedSequence">IAlignedSequence instance.</param>
        internal PairwiseAlignedSequence(IAlignedSequence alignedSequence)
            : base(alignedSequence)
        {
            // Impl.
        }
        #endregion

        #region Properties
        /// <summary>
        /// Gets or sets Alignment of First Sequence.
        /// </summary>
        public ISequence FirstSequence
        {
            get
            {
                if (Sequences.Count > 0)
                {
                    return Sequences[0];
                }

                return null;
            }

            set
            {
                if (Sequences.Count == 0)
                {
                    Sequences.Add(null);
                }

                Sequences[0] = value;
            }
        }

        /// <summary>
        /// Gets or sets Alignment of Second Sequence.
        /// </summary>
        public ISequence SecondSequence
        {
            get
            {
                if (Sequences.Count > 1)
                {
                    return Sequences[1];
                }

                return null;
            }

            set
            {
                if (Sequences.Count == 0)
                {
                    Sequences.Add(null);
                }

                if (Sequences.Count == 1)
                {
                    Sequences.Add(null);
                }

                Sequences[1] = value;
            }
        }

        /// <summary>
        /// Gets or sets Consensus of FirstSequence and SecondSequence.
        /// </summary>
        public ISequence Consensus
        {
            get
            {
                if (Metadata.ContainsKey(ConsensusKey))
                {
                    return Metadata[ConsensusKey] as ISequence;
                }

                return null;
            }

            set
            {
                Metadata[ConsensusKey] = value;
            }
        }

        /// <summary>
        /// Gets or sets Score of the alignment.
        /// </summary>
        public long Score
        {
            get
            {
                long score = 0;
                if (Metadata.ContainsKey(ScoreKey))
                {
                    if (Metadata[ScoreKey] is long)
                    {
                        score = (long)Metadata[ScoreKey];
                    }
                }

                return score;
            }

            set
            {
                Metadata[ScoreKey] = value;
            }
        }

        /// <summary>
        /// Gets or sets Offset of FirstSequence.
        /// 
        /// Danger!  This is not the start position of the alignment, it is some goofy
        /// difference in start positions.
        /// </summary>
        public long FirstOffset
        {
            get
            {
                long score = 0;
                if (Metadata.ContainsKey(FirstOffsetKey))
                {
                    if (Metadata[FirstOffsetKey] is long)
                    {
                        score = (long)Metadata[FirstOffsetKey];
                    }
                }

                return score;
            }

            set
            {
                Metadata[FirstOffsetKey] = value;
            }
        }
        /// <summary>
        /// Converts a query positiion (0-indexed) to a refernce position, or returns null if no overlap.
        /// For example:
        ///         0123456     7
        /// Ref:    AAAAAAA-----AAA----
        /// Query:  AAAAAAACCCCCA
        ///         0123456789012
        /// 
        /// Returns 12 when given .  Note: it will account for clipping if this is a local alignment (that is it will return global coordinates).
        /// </summary>
        /// <returns>The query position correspondingto reference position.</returns>
        /// <param name="refPos">Reference position.</param>
        public long? FindQueryPositionCorrespondingtoReferencePosition(int refPos)
        {
            if (!FirstSequenceStart.HasValue || !SecondSequenceStart.HasValue) {
                throw new BioinformaticsException ("");
            }

            var gap = DnaAlphabet.Instance.Gap;
            int r_pos = (int) FirstSequenceStart.Value - 1;
            int q_pos = (int)SecondSequenceStart.Value - 1;

            int start = r_pos;
            int end = r_pos + (int)FirstSequence.Count;
            if (refPos < start || refPos > end) {
                return null;
            }

            for (int i = 0; i < FirstSequence.Count; i++) {
                if (FirstSequence [i] != gap) {
                    r_pos++;
                }
                if (SecondSequence [i] != gap) {
                    q_pos++;
                }
                if (r_pos == refPos) {
                    if (SecondSequence [i] != gap) {
                        return q_pos;
                    }
                    break;
                }
            }
            return null;       
        }

      

        /// <summary>
        /// The first 0 indexed base of the alignment.  INCLUSIVE
        /// </summary>
        public long? FirstSequenceStart;
        public long? SecondSequenceStart;

        /// <summary>
        /// Gets or sets Offset of SecondSequence.
        /// Danger!  This is not the start position of the alignment, it is some goofy
        /// difference in start positions.
        /// </summary>
        public long SecondOffset
        {
            get
            {
                long score = 0;
                if (Metadata.ContainsKey(SecondOffsetKey))
                {
                    if (Metadata[SecondOffsetKey] is long)
                    {
                        score = (long)Metadata[SecondOffsetKey];
                    }
                }

                return score;
            }

            set
            {
                Metadata[SecondOffsetKey] = value;
            }
        }
        #endregion
        /// <summary>
        /// Converts the Consensus, First and Second sequences.
        /// </summary>
        /// <returns>Consensus, First and Second sequences.</returns>
        public override string ToString()
        {
            StringBuilder builder = new StringBuilder();
            builder.AppendLine ("Alignment to " + this.FirstSequence.ID);
            builder.AppendLine ("Start Ref = " + this.FirstSequenceStart.ToString () + " ; Start Query = " + this.SecondSequenceStart.ToString ());
            if (Consensus != null) {
                builder.AppendLine (this.Consensus.ConvertToString ());
            }
            builder.AppendLine(this.FirstSequence.ConvertToString());
            builder.AppendLine(this.SecondSequence.ConvertToString());
            return builder.ToString();
        }
    }
}
