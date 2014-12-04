using System;
using System.Linq;
using System.Collections.Generic;
using PacBio.IO;
using PacBio.Utils;

namespace PacBio.Consensus
{
    public enum MutationType {INSERTION, DELETION, SUBSTITUTION}
    /// <summary>
    /// Used to represent a potential modification to the consensus template.  Typically a Mutation always applies to
    /// the forward strand.  The Strand convert method converts the mutation to apply to the reverse strand this usage
    /// should be kept as localized as possible.
    /// </summary>
    public struct Mutation : IComparable<Mutation>
    {
        public int TemplatePosition;
        public MutationType Type;
        public char Base;

        // Apply a mutation to a template
        public string Apply(string template)
        {
            var tt = template.ToCharArray().ToList();

            switch (Type)
            {
                case MutationType.INSERTION:
                    tt.Insert(TemplatePosition, Base);
                    return new string(tt.ToArray());

                case MutationType.DELETION:
                    tt.RemoveAt(TemplatePosition);
                    return new string(tt.ToArray());

                case MutationType.SUBSTITUTION:
                    tt[TemplatePosition] = Base;
                    return new string(tt.ToArray());

                default:
                    throw new ArgumentException("Unrecognized mutation type");
            }
        }

        /// <summary>
        /// Apply a set of mutations to a template. Mutations must be separated by at least two 
        /// base, or weird interactions may occur
        /// </summary>
        /// <param name="mutations">Set of mutations to apply</param>
        /// <param name="template">TrialTemplate to mutate</param>
        public static string ApplyMany(List<Mutation> mutations, string template)
        {
            mutations.Sort();
            mutations.Reverse();

            var tt = template.ToCharArray().ToList();

            foreach (Mutation m in mutations)
            {
                switch (m.Type)
                {
                    case MutationType.INSERTION:
                        tt.Insert(m.TemplatePosition, m.Base);
                        break;

                    case MutationType.DELETION:
                        tt.RemoveAt(m.TemplatePosition);
                        break;

                    case MutationType.SUBSTITUTION:
                        tt[m.TemplatePosition] = m.Base;
                        break;

                    default:
                        throw new ArgumentException("Unrecognized mutation type");
                }
            }
            return new string(tt.ToArray());
        }

        public int CompareTo(Mutation other)
        {
            return TemplatePosition.CompareTo(other.TemplatePosition);
        }

        /// <summary>
        /// Check whether this mutation will leave the template unchanged.
        /// </summary>
        /// <param name="template">The template that will be mutated, in channel space</param>
        public bool IsSynonymous(string template)
        {
            // Insertions and deletions are never synonymous,
            // mismatches are synonymous if the target channel matches the current channel
            switch (Type)
            {
                case MutationType.INSERTION:
                    return false;

                case MutationType.DELETION:
                    return false;

                case MutationType.SUBSTITUTION:
                    return template[TemplatePosition] == Base;

                default:
                    throw new ArgumentException("Unrecognized mutation type");
            }
        }

        public bool IsSynonymous(TrialTemplate trialTemplate)
        {
            return IsSynonymous(trialTemplate.GetSequence(Strand.Forward));
        }

        /// <summary>
        /// Convert the mutation into the matching mutation on another strand
        /// </summary>
        /// <param name="strand">The strand to conver the mutation to</param>
        /// <param name="trialTemplate">The template to which the mutation will be applied</param>
        /// <returns>A mutation object converted to the target strand</returns>
        public Mutation StrandConvert(Strand strand, TrialTemplate trialTemplate)
        {
            if (strand == Strand.Forward)
                return this;

            var newBase = Base.Complement();
            int newPos;

            if (Type == MutationType.INSERTION)
            {
                newPos = trialTemplate.Sequence.Length - TemplatePosition;
            }
            else
            {
                newPos = trialTemplate.Sequence.Length - TemplatePosition - 1;
            }

            return new Mutation { Base = newBase, TemplatePosition = newPos, Type = Type };
        }

        public override string ToString()
        {
            return String.Format("Pos: {0}, Mode: {1}, Base: {2}", TemplatePosition, Type, Base);
        }

        public int LengthDiff
        {
            get
            {
                switch (Type)
                {
                    case MutationType.SUBSTITUTION:
                        return 0;
                    case MutationType.DELETION:
                        return -1;
                    case MutationType.INSERTION:
                        return 1;
                    default:
                        throw new Exception("unrecognized MutationType");
                }
            }
        }

    }
    
    /// <summary>
    /// Whether to recurse using the forward-backward algorithm or the viterbi algorithm
    /// </summary>
    public enum RecursionAlgo
    {
        Viterbi=0,
        Prob=1
    }
}
