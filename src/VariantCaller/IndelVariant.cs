using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Bio;

namespace VariantCaller
{
    public enum IndelType : byte { Insertion, Deletion}
    public class IndelVariant : Variant
    {
        /// <summary>
        /// Insertion or Deletion
        /// </summary>
        public readonly IndelType InsertionOrDeletion;

        /// <summary>
        /// These are the bases in the insertion (or what was deleted).
        /// </summary>
        public readonly string InsertedOrDeletedBases;
        public IndelVariant(int position, int length, Sequence reference, string bases, IndelType insertionOrDeletion, bool atAlignmentEnd = false) 
            : base (position, reference, atAlignmentEnd)
        {
            Type = VariantType.INDEL;
            InsertionOrDeletion = insertionOrDeletion;
            Length = length;
            InsertedOrDeletedBases = bases;
        }
    }
}
