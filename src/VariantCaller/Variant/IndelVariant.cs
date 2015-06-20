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


        /// <summary>
        /// True if indel has back to back bases.
        /// </summary>
        /// <value><c>true</c> if in homo polymer; otherwise, <c>false</c>.</value>
        public bool InHomopolymer {
            get {
                if (_homopolymerLength == -1) {
                    _homopolymerLength = determineHomoPolymerLength ();
                }
                return _homopolymerLength > 1;
            }
        }

        private int _homopolymerLength = -1;
        ///
        /// Number of bases after the
        ///
        public int HomopolymerLengthInReference {
            get{
                if (_homopolymerLength == -1) {
                    _homopolymerLength = determineHomoPolymerLength ();
                }
                return _homopolymerLength;
            }

        }

        private char _homopolymerBase;

        public char HomopolymerBase
        {
            get {
                if (_homopolymerLength == -1) {
                    _homopolymerLength = determineHomoPolymerLength ();
                }
                return _homopolymerBase;
            }
        }

        public IndelVariant(int position, int length, ISequence reference, string bases, IndelType insertionOrDeletion, bool atAlignmentEnd = false) 
            : base (position, reference, atAlignmentEnd)
        {
            Type = VariantType.INDEL;
            InsertionOrDeletion = insertionOrDeletion;
            Length = length;
            InsertedOrDeletedBases = bases;
        }

        /// <summary>
        /// Returns the length of the homo polymer, 1 if 
        /// </summary>
        /// <returns>The homo polymer length.</returns>
        private int determineHomoPolymerLength()
        {
            if (this.AtEndOfAlignment) {
                return int.MaxValue;
            }

            int pos = this.StartPosition + 1;
            // Since this is left aligned, anything onwards is considered a homopolymer
            if (pos > RefSeq.Count) {
                throw new BioinformaticsException (" Somehow ran off the reference here, error");
            }
            byte start_bp = this.RefSeq [pos];
            _homopolymerBase = (char)start_bp;
            int len = 0;
            while (pos < RefSeq.Count && RefSeq[pos] == start_bp) {
                len++;
                pos++;
            }
            return len;
        }
        public override string ToString ()
        {
            var insert = this.InsertionOrDeletion == IndelType.Deletion ? "Deletion" : "Insertion";
            return string.Format ("[IndelVariant: InHomopolymer={0}, HomopolymerLengthInReference={1}, HomopolymerBase={2}, Type={3}]", InHomopolymer, HomopolymerLengthInReference, HomopolymerBase, insert);
        }
    }
}
