using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Bio;
using Bio.Algorithms.MUMmer;
using Bio.Algorithms.Alignment;
using Bio.Extensions;

namespace VariantCaller
{
    public class NucmerQueryable
    {
        private NUCmer nucmer;
        private Sequence ReferenceSequence;
        #region NUCMER STUFF

        /// <summary>
        /// Use anchor matches that are unique in both the reference and query.
        /// </summary>
        public static bool Mum = false;

        /// <summary>
        /// Use anchor matches that are unique in the reference but not necessarily unique in the query (default behavior).
        /// </summary>
        public static bool MumReference = true;

        /// <summary>
        /// Use all anchor matches regardless of their uniqueness.
        /// </summary>
        public static bool MaxMatch = true;

        /// <summary>
        /// Distance an alignment extension will attempt to extend poor scoring regions before giving up (default 200).
        /// </summary>
        public static int BreakLength = 200;

        /// <summary>
        /// Minimum cluster length (default 65).
        /// </summary>
        public static int MinCluster = 35;
        /// <summary>
        /// Maximum diagonal difference factor for clustering, i.e. diagonal difference / match separation (default 0.12).
        /// </summary>
        public static double DiagFactor = 0.12;

        /// <summary>
        /// Maximum gap between two adjacent matches in a cluster (default 90).
        /// </summary>
        public static int MaxGap = 90;

        /// <summary>
        /// Minimum length of an maximal exact match (default 20).
        /// </summary>
        public static int MinMatch = 20;

        /// <summary>
        /// Align only the reverse strand of the query sequence to the forward strand of the reference.
        /// </summary>
        public static bool Reverse = false;

        /// <summary>
        /// Toggle the outward extension of alignments from their anchoring clusters. Setting --noextend will 
        /// prevent alignment extensions but still align the DNA between clustered matches and create the .delta file 
        /// (default --extend).
        /// </summary>
        public static bool NotExtend = false;

        /// <summary>
        ///  Gets or sets maximum fixed diagonal difference
        /// </summary>
        public static int FixedSeparation = 5;
        #endregion 
        public NucmerQueryable(Sequence touse)
        {
            ReferenceSequence = touse;
            nucmer= new NUCmer(ReferenceSequence)
          {
                    FixedSeparation = FixedSeparation,
                    BreakLength = BreakLength,
                    LengthOfMUM = MinMatch,
                    MaximumSeparation = MaxGap,
                    MinimumScore = MinCluster,
                    SeparationFactor = (float) DiagFactor
                };
        }
        public  IEnumerable<IEnumerable<DeltaAlignment>> GetDeltaAlignments(ISequence querySequence)
        {
            IEnumerable<ISequence> querySequences = AddReverseComplementsToSequenceList(querySequence);
            foreach (ISequence qs in querySequences)
            {
                yield return nucmer.GetDeltaAlignments(qs, !MaxMatch, qs.IsMarkedAsReverseComplement());
            }
            
        }
        /// <summary>
        /// Given a list of sequences, create a new list with the original sequence followed
        /// by the Reverse Complement of that sequence.
        /// </summary>
        /// <param name="sequenceList">List of sequence.</param>
        /// <returns>Returns the List of sequence.</returns>
        private IEnumerable<ISequence> AddReverseComplementsToSequenceList(ISequence seq)
        {
            yield return seq;

            ISequence rcSequence = seq.GetReverseComplementedSequence();
            if (rcSequence != null)
            {
                rcSequence.MarkAsReverseComplement();
                yield return rcSequence;
            }

        }
    }
}
