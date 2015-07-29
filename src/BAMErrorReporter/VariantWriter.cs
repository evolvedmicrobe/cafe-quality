using System;
using Bio;
using Bio.IO.BAM;
using Bio.IO.SAM;
using System.Collections.Generic;
using System.Linq;
using Bio.IO.PacBio;
using System.Diagnostics;
using System.IO;
using VariantCaller;

namespace BAMErrorReporter
{
    public class VariantWriter
    {
        public StreamWriter sw;
        public VariantWriter (string filename)
        {
            sw = new StreamWriter (filename);
            sw.WriteLine ("Ref,Pos,zmw,type,length,homopolymerLength,homopolymerChar,indelSize,indeltype,QV");

        }

        public void Write(PacBioCCSRead read, Variant variant) {
            var hplength = homopolymerLength (variant);
            var vtype = variant.Type.ToString ();
            string indeltype = "NA";
            string indelLength = "NA";
            string homoChar = "N";
            if (!variant.AtEndOfAlignment && variant.Type == VariantType.INDEL) {
                var vi = (variant as IndelVariant);
                indeltype = vi.InsertionOrDeletion.ToString ();
                indelLength = vi.InsertedOrDeletedBases.Length.ToString();
                homoChar = vi.HomopolymerBase.ToString ();
            }
            var toW = String.Join (",", variant.RefSeq.ID,
                                        variant.StartPosition,
                                        read.HoleNumber,
                                        vtype,
                                        read.Sequence.Count,
                                        hplength,
                                        homoChar,
                                        indelLength,
                                        indeltype,
                                        variant.QV.ToString ());
            sw.WriteLine (toW);
            sw.Flush ();
        }

        private string homopolymerLength(Variant v)
        {
            if (v.AtEndOfAlignment)
                return "-999";
            if (v is IndelVariant) {
                return (v as IndelVariant).HomopolymerLengthInReference.ToString ();
            } else if (v is SNPVariant) {
                return "1";
            }
            throw new Exception ("Unknown variant type");

        }

    }
}

