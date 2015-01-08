using System;
using System.Text;
using System.Linq;
using Bio;
using PacBio.Data;
using PacBio.Utils;
using System.Diagnostics;
using System.Collections.Generic;

namespace VariantCaller
{
    /// <summary>
    /// This is a ridiculous class that only exists until I understand the data model better.
    /// It appears that the BaxH5Reader produces reads which can produce data when requested from
    /// a method.  This gets around this by 
    /// </summary>
    public class ReadFromZMW
    {

        [OutputAttribute]
        public string AssignedReferenceName;

        private Reference assignedReference;

        public Reference AssignedReference
        {
            get{ return assignedReference; }
            set { assignedReference = value;
                AssignedReferenceName = value.RefSeq.ID;
            }
        }

        [OutputArrayAttribute]
        public readonly byte[] DeletionQV;
        [OutputArrayAttribute]
        public readonly byte[] DeletionTag;
        [OutputArrayAttribute]
        public readonly short[] IpdInFrames;
        [OutputAttribute]
        public readonly string BaseCalls;
        [OutputArrayAttribute]
        public readonly byte[] InsertionQV;
        [OutputArrayAttribute]
        public readonly byte[] MergeQV;
        [OutputArrayAttribute]
        public readonly byte[] SubstitutionQV;
        [OutputArrayAttribute]
        public readonly short[] PulseWidthInFrames;

        //TODO: These should all be removed later.
        [OutputAttribute]
        public bool ReverseComplementedOriginally;
        [OutputAttribute]
        public string ConsensusIndelSize;
        [OutputAttribute]
        public float NoErrorSumProductScoreP6;
        [OutputAttribute]
        public float NoErrorViterbiScoreP6;
        [OutputAttribute]
        public float OneDeletionSumProductScoreP6;
        [OutputAttribute]
        public float OneDeletionErrorViterbiScoreP6;

        [OutputAttribute]
        public float NoErrorSumProductScoreC2;
        [OutputAttribute]
        public float NoErrorViterbiScoreC2;
        [OutputAttribute]
        public float OneDeletionSumProductScoreC2;
        [OutputAttribute]
        public float OneDeletionErrorViterbiScoreC2;


        [OutputAttribute]
        public float FullLengthCorrectSumProductScore;
        [OutputAttribute]
        public float FullLengthIncorrectSumProductScore;
        [OutputAttribute]
        public int OriginalSubReadLength;
        [OutputAttribute]
        public int SubReadNumber;
        [OutputAttribute]
        public int HPSectionLength;

        [OutputAttribute]
        public float RQ;
        [OutputAttribute]
        public int Zmw;
        /// <summary>
        /// IF we do an alignment, how long this is.
        /// </summary>
        [OutputAttribute]
        public int AlignedLength;

        /// <summary>
        /// The total count of DelTags not equal to 'N'
        /// </summary>
        [OutputAttribute]
        public int CountDelTags;
        /// <summary>
        /// The count of DelTags that could be considered correct
        /// </summary>
        [OutputAttribute]
        public int CountCorrectDelTags;
        [OutputAttribute]
        public int SpikeMergeQVCount;
        [OutputAttribute]
        public int HPDeletionTagNotNCount;
        [OutputAttribute]
        public string HomopolymerLengthFromAlignment;
        [OutputAttribute]
        public int HomopolymerLengthFromCounting;
        [OutputAttribute]
        public int SubRead;
        [OutputAttribute]
        public float SummedPulseWidthForHP;


        public ReadFromZMW()
        {
        }
        // End output section

        public ReadFromZMW (ZmwRead baseRead )
        {

            var features = baseRead.Features;
            DeletionQV = features.DeletionQV ();
            IpdInFrames = features.IpdInFrames();
            BaseCalls = features.Basecalls ();
            InsertionQV = features.InsertionQV();
            MergeQV = features.MergeQV ();
            SubstitutionQV = features.SubstitutionQV ();
            PulseWidthInFrames = features.PulseWidthInFrames ();
            DeletionTag = features.DeletionTag ();
        }

        ReadFromZMW(byte[] deletionQV, short[] ipdInFrames, string baseCalls, 
            byte[] insertionQV, byte[] mergeQV, byte[] substitutionQV, short[] pulseWidthInFrames, byte[] deletionTag)
        {
            this.DeletionQV = deletionQV;
            this.IpdInFrames = ipdInFrames;
            this.BaseCalls = baseCalls;
            this.InsertionQV = insertionQV;
            this.MergeQV = mergeQV;
            this.SubstitutionQV = substitutionQV;
            this.PulseWidthInFrames = pulseWidthInFrames;
            this.DeletionTag = deletionTag;
        }


        public ReadFromZMW GetSubSection(int start, int length)
        {
            // Validate match
            if (length + start > this.BaseCalls.Length) {
                throw new ArgumentOutOfRangeException ();
            }
            // Slice out relevant bits

            var dqv = DeletionQV.Slice (start, length);
            var mqv = MergeQV.Slice (start, length);
            var ipdif = IpdInFrames.Slice (start, length);
            var bc = BaseCalls.Substring (start, length);
            var iqv = InsertionQV.Slice (start, length);
            var subqv = SubstitutionQV.Slice (start, length);
            var pwif = PulseWidthInFrames.Slice (start, length);
            var delt = DeletionTag.Slice (start, length);

            // Return new subread
            return new ReadFromZMW (dqv, ipdif, bc, iqv, mqv, subqv, pwif, delt);
        }

        public ReadFromZMW GetSubRead(CCSSubRead sub)
        {
            // Validate match
            int len = (int)sub.End - (int)sub.Start;
            Trace.Assert (len == sub.Seq.Length);
            var seq = new char[len];
            for (int i = 0; i < seq.Length; i++) {
                var cur = (char) sub.Seq [i];
                var other = this.BaseCalls [i + (int)sub.Start];
                Trace.Assert (cur == other);
                seq [i] = cur;
            }
            var start = (int)sub.Start;
            return GetSubSection (start, len);
        }

        public static List<string> CSVHeaderFields()
        {
            return OutputHelper.GetHeaders(typeof(ReadFromZMW));
        }
        public ReadFromZMW GetReverseComplement()
        {
            var s = new Sequence (DnaAlphabet.Instance, BaseCalls, false);
            var ds = new Sequence (AmbiguousDnaAlphabet.Instance, DeletionTag, false);
            var ds_rc = ds.GetReverseComplementedSequence ().ConvertToString ().Select (z => (byte)z).ToArray ();
            return new ReadFromZMW (DeletionQV.Reverse ().ToArray (),
                IpdInFrames.Reverse ().ToArray (),
                s.GetReverseComplementedSequence ().ConvertToString (),
                InsertionQV.Reverse ().ToArray (),
                MergeQV.Reverse ().ToArray (),
                SubstitutionQV.Reverse ().ToArray (),
                PulseWidthInFrames.Reverse ().ToArray (),
                ds_rc);                         
        }
         
    }
}

