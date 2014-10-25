using System;
using System.Linq;
using PacBio.Data;

namespace VariantCaller
{
    /// <summary>
    /// This is a ridiculous class that only exists until I understand the data model better.
    /// It appears that the BaxH5Reader produces reads which can produce data when requested from
    /// a method.  This gets around this by 
    /// </summary>
    public class ReadFromZMW
    {
        [OutputArrayAttribute]
        public readonly byte[] DeletionQV;
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
        }

        public static void CSVHeaderFields()
        {
          
        }
    }
}

