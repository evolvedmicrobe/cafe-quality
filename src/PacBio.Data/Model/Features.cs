using System;
using System.Diagnostics;
using PacBio.Utils;

namespace PacBio.Data
{
    /// <summary>
    /// An ISequenceFeatures is a "data frame" representing a sequence of PacBio basecalls
    /// and the additional rich pulse features (QVs and kinetic data).
    /// The order of the bases/features is always the order the bases were read by
    /// the sequencing instrument ("native" order).
    /// 
    /// Implementing classes come in two flavors:
    ///   - "Cursor" types, which are tied to a data file and generally may only
    ///     be used from a single thread
    ///   - "Manifest" (fully-in-memory) types, which can be used freely and should
    ///     be used as the data type for algorithms to enable parallelization.  Manifest
    ///     types, defined in this file, implement a copy constructor that enables creation
    ///     from corresponding cursor types.
    /// </summary>
    public interface ISequenceFeatures
    {
        int Length { get; }
        string Basecalls();
    }

    public interface IQvFeatures : ISequenceFeatures
    {
        byte[] InsertionQV();
        byte[] DeletionQV();
        byte[] SubstitutionQV();
        byte[] MergeQV();
        byte[] DeletionTag();
    }

    public interface IKineticsFeatures : ISequenceFeatures
    {
        short[] IpdInFrames();
        short[] PulseWidthInFrames();
    }

    public interface IFullFeatures : IQvFeatures, IKineticsFeatures {}

    public class SequenceFeatures : ISequenceFeatures
    {
        protected string basecalls;

        public int Length { get { return basecalls.Length; } }
        public string Basecalls() { return basecalls; }

        public SequenceFeatures(string basecalls)
        {
            this.basecalls = basecalls;
        }

        public SequenceFeatures(ISequenceFeatures other)
            : this(other.Basecalls())
        {}

        protected SequenceFeatures() {}

        public virtual SequenceFeatures Slice(int b, int l)
        {
            return new SequenceFeatures { basecalls = basecalls.Substring(b, l) };
        }
    }

    public class QvFeatures : SequenceFeatures, IQvFeatures
    {
        protected byte[] insertionQV;
        protected byte[] deletionQV;
        protected byte[] substitutionQV;
        protected byte[] mergeQV;
        protected byte[] deletionTag;

        public QvFeatures(IQvFeatures other)
            :base(other)
        {
            insertionQV    = other.InsertionQV();
            deletionQV     = other.DeletionQV();
            substitutionQV = other.SubstitutionQV();
            mergeQV        = other.MergeQV();
            deletionTag    = other.DeletionTag();
        }

        public byte[] InsertionQV()    { return insertionQV;    }
        public byte[] DeletionQV()     { return deletionQV;     }
        public byte[] SubstitutionQV() { return substitutionQV; }
        public byte[] MergeQV()        { return mergeQV;        }
        public byte[] DeletionTag()    { return deletionTag;    }

        protected QvFeatures() {}

        public override SequenceFeatures Slice(int b, int l)
        {
            return new QvFeatures() {
                basecalls      = basecalls      .Substring(b, l),
                insertionQV    = insertionQV    .Slice(b, l),
                deletionQV     = deletionQV     .Slice(b, l),
                substitutionQV = substitutionQV .Slice(b, l),
                mergeQV        = mergeQV        .Slice(b, l),
                deletionTag    = deletionTag    .Slice(b, l)
            };
        }
    }

    public class KineticsFeatures : SequenceFeatures, IKineticsFeatures
    {
        short[] ipdInFrames;
        short[] pulseWidthInFrames;

        public KineticsFeatures(IKineticsFeatures other)
            :base(other)
        {
            ipdInFrames = other.IpdInFrames();
            pulseWidthInFrames = other.PulseWidthInFrames();
        }

        protected KineticsFeatures() {}

        public override SequenceFeatures Slice(int b, int l)
        {
            return new KineticsFeatures() {
                basecalls          = basecalls          .Substring(b, l),
                ipdInFrames        = ipdInFrames        .Slice(b, l),
                pulseWidthInFrames = pulseWidthInFrames .Slice(b, l)
            };
        }

        public short[] IpdInFrames()        { return ipdInFrames; }
        public short[] PulseWidthInFrames() { return pulseWidthInFrames; }
    }


    public class FullFeatures : SequenceFeatures, IFullFeatures
    {
        byte[] insertionQV;
        byte[] deletionQV;
        byte[] substitutionQV;
        byte[] mergeQV;
        byte[] deletionTag;
        short[] ipdInFrames;
        short[] pulseWidthInFrames;

        public FullFeatures(IFullFeatures other)
            :base(other)
        {
            insertionQV    = other.InsertionQV();
            deletionQV     = other.DeletionQV();
            substitutionQV = other.SubstitutionQV();
            mergeQV        = other.MergeQV();
            deletionTag    = other.DeletionTag();
            ipdInFrames = other.IpdInFrames();
            pulseWidthInFrames = other.PulseWidthInFrames();
        }

        protected FullFeatures() {}

        public override SequenceFeatures Slice(int b, int l)
        {
            return new FullFeatures() {
                basecalls          = basecalls          .Substring(b, l),
                insertionQV        = insertionQV        .Slice(b, l),
                deletionQV         = deletionQV         .Slice(b, l),
                substitutionQV     = substitutionQV     .Slice(b, l),
                mergeQV            = mergeQV            .Slice(b, l),
                deletionTag        = deletionTag        .Slice(b, l),
                ipdInFrames        = ipdInFrames        .Slice(b, l),
                pulseWidthInFrames = pulseWidthInFrames .Slice(b, l)
            };
        }

        public byte[] InsertionQV()    { return insertionQV;    }
        public byte[] DeletionQV()     { return deletionQV;     }
        public byte[] SubstitutionQV() { return substitutionQV; }
        public byte[] MergeQV()        { return mergeQV;        }
        public byte[] DeletionTag()    { return deletionTag;    }
        public short[] IpdInFrames()        { return ipdInFrames; }
        public short[] PulseWidthInFrames() { return pulseWidthInFrames; }
    }

}
