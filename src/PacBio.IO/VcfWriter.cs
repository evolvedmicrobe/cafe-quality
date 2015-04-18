using System;
using System.Collections.Generic;
using System.Collections.Specialized;
using System.IO;
using System.Linq;
using System.Security.Cryptography;
using System.Text;
using PacBio.Utils;

namespace PacBio.IO
{
    public struct Variant
    {
        /*
        public string refName;
        public int    position;
        public string refValue;
        public string altValue;
         */
        public int TemplatePosition;
        public char RefBase;
        public char AltBase;
        public double Fraction;
        public double FractionLow;
        public double FractionHigh;
        public double PValue;
        public double PseudoCounts;
        public double Coverage;
    }


    public class VcfContig
    {
        public readonly string Id;
        public readonly int Length;
        public readonly string Md5;
        public readonly string Sequence;
        public readonly List<Variant> Variants;
        public readonly bool AminoVariants;

        public VcfContig(string id, string sequence, bool aminoVariants = false)
        {
            Id = id;
            Length = sequence.Length;
            Md5 = ComputeStringMd5(sequence);
            Sequence = sequence;
            Variants = new List<Variant>();
            AminoVariants = aminoVariants;
        }

        public void AddVariant(Variant variant)
        {
            Variants.Add(variant);
        }

        public void AddVariants(IEnumerable<Variant> variants)
        {
            Variants.AddRange(variants);
        }

        private static string ComputeStringMd5(string value)
        {
            using (var md5 = MD5.Create())
            {
                return BitConverter
                    .ToString(md5.ComputeHash(Encoding.ASCII.GetBytes(value)))
                    .Replace("-", string.Empty)
                    .ToLower();
            }
        }
    }


    abstract public class VariantWriter : IDisposable
    {
        protected readonly FileStream _fileStream;
        protected readonly StreamWriter _streamWriter;
        private readonly OrderedDictionary _contigs;
        protected readonly string _referenceFilename;

        protected VariantWriter(string filename, string referenceFilename)
        {
            _fileStream = new FileStream(filename, FileMode.Create);
            _streamWriter = new StreamWriter(_fileStream);
            _contigs = new OrderedDictionary();
            _referenceFilename = referenceFilename;
        }

        public VcfContig AddContig(string id, string sequence, bool aminoVariants = false)
        {
            var contig = new VcfContig(id, sequence, aminoVariants);

            _contigs.Add(id, contig);

            return contig;
        }

        protected abstract void WriteContents(ICollection<VcfContig> contigs);

        public void Dispose()
        {
            Dispose(true);
            GC.SuppressFinalize(this);
        }

        ~VariantWriter()
        {
            Dispose(false);
        }

        private bool disposed = false;

        protected virtual void Dispose(bool disposing)
        {
            WriteContents(_contigs.Values.OfType<VcfContig>().ToList());

            if (disposed)
                return;

            _streamWriter.Close();
            _fileStream.Close();

            if (disposing)
            {
                _streamWriter.Dispose();
                _fileStream.Dispose();
            }

            disposed = true;
        }

        protected static string AminoVariant(string dna, int pos, char altBase)
        {
            var tplStart = 3 * (pos / 3);

            var codon = Enumerable.Range(tplStart, 3)
                .Select(p => dna[p])
                .ToArray();

            var refCodon = String.Concat(codon);

            codon[pos % 3] = altBase;
            var altCodon = String.Concat(codon);

            return String.Format("{0}{1}{2}",
                                 Translation.TranslateCodon(refCodon),
                                 pos / 3 + 1,
                                 Translation.TranslateCodon(altCodon));
        }

        protected static int ErrorProbabilityToQv(double value)
        {
            return value == 0 ? 93 : Math.Min(93, (int)Math.Round(-10 * Math.Log10(value)));
        }
    }


    public class CsvVariantWriter : VariantWriter
    {
        private string _delimiter;

        public CsvVariantWriter(string filename, string referenceFilename, string delimiter = ",")
            : base(filename, referenceFilename)
        {
            _delimiter = delimiter;
        }

        protected override void WriteContents(ICollection<VcfContig> contigs)
        {
            var aminoVariants = contigs.Aggregate(false, (av, contig) => av || contig.AminoVariants);

            var header = new List<String>
                         {
                             "CONTIG", "POSITION", "REF", "ALT", "COVERAGE", "SUPPORT",
                             "FREQUENCY", "\"FREQUENCY LOW\"", "\"FREQUENCY HIGH\"", "\"P-VALUE\""
                         };

            if (aminoVariants)
            {
                header.Add("\"AMINO ACID TRANSLATION\"");
            }

            _streamWriter.WriteLine(string.Join(_delimiter, header));

            foreach (var contig in contigs)
            {
                foreach (var variant in contig.Variants)
                {
                    var line = new List<String>
                               {
                                   String.Format("\"{0}\"", contig.Id),
                                   String.Format("{0}", variant.TemplatePosition + 1),
                                   String.Format("{0}", variant.RefBase),
                                   String.Format("{0}", variant.AltBase),
                                   String.Format("{0}", (int) Math.Floor(variant.Coverage)),
                                   String.Format("{0}", variant.PseudoCounts),
                                   String.Format("{0:0.000}", variant.Fraction),
                                   String.Format("{0:0.000}", variant.FractionLow),
                                   String.Format("{0:0.000}", variant.FractionHigh),
                                   String.Format("{0}", ErrorProbabilityToQv(variant.PValue))
                               };

                    if (contig.AminoVariants)
                    {
                        line.Add(AminoVariant(contig.Sequence,
                                              variant.TemplatePosition,
                                              variant.AltBase));
                    }

                    _streamWriter.WriteLine(string.Join(_delimiter, line));
                }
            }
        }
    }


    public class VcfWriter : VariantWriter
    {
        public VcfWriter(string filename, string referenceFilename)
            : base(filename, referenceFilename)
        { }

        protected override void WriteContents(ICollection<VcfContig> contigs)
        {
            // VCF spec requires UNIX line endings
            _streamWriter.NewLine = "\n";

            // write the header
            _streamWriter.WriteLine("##fileformat=VCFv4.1");
            _streamWriter.WriteLine("##source=ConsensusTools");
            _streamWriter.WriteLine(String.Format("##reference={0}", _referenceFilename));

            var aminoVariants = false;

            foreach (var contig in contigs)
            {
                _streamWriter.WriteLine(String.Format("##contig=<ID={0},length={1},md5={2}>",
                        contig.Id, contig.Length, contig.Md5));
                aminoVariants |= contig.AminoVariants;
            }

            _streamWriter.WriteLine("##INFO=<ID=CV,Number=1,Type=Integer,Description=\"Coverage\">");
            _streamWriter.WriteLine("##INFO=<ID=ZS,Number=A,Type=Integer,Description=\"ZMW support\">");
            _streamWriter.WriteLine("##INFO=<ID=VF,Number=A,Type=Float,Description=\"Variant frequency\">");
            _streamWriter.WriteLine("##INFO=<ID=CL,Number=A,Type=Float,Description=\"Lower confidence interval\">");
            _streamWriter.WriteLine("##INFO=<ID=CH,Number=A,Type=Float,Description=\"Upper confidence interval\">");
            _streamWriter.WriteLine("##INFO=<ID=PV,Number=A,Type=Integer,Description=\"P-value (Phred-encoded)\">");

            if (aminoVariants)
            {
                _streamWriter.WriteLine("##INFO=<ID=AA,Number=A,Type=String,Description=\"Amino acid translation\">");
            }

            _streamWriter.WriteLine("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");

            // write the variants themselves
            foreach (var contig in contigs)
            {
                WriteVariantsForContig(contig);
            }
        }

        private void WriteVariantsForContig(VcfContig contig)
        {
            var variantsBySite = contig.Variants
                .GroupBy(variant => variant.TemplatePosition)
                .Select(variantGroup => variantGroup.ToList());

            foreach (var variantGroup in variantsBySite)
            {
                if (!variantGroup.Any())
                    continue;

                var fst = variantGroup.First();
                var pos = fst.TemplatePosition;
                var refBase = fst.RefBase;
                var altBases = String.Join(",", variantGroup.Select(mv => mv.AltBase));
                var cov = fst.Coverage;
                var counts = String.Join(",", variantGroup.Select(mv => mv.PseudoCounts));
                var fracs = String.Join(",", variantGroup.Select(mv => mv.Fraction.ToString("0.000")));
                var fracLows = String.Join(",", variantGroup.Select(mv => mv.FractionLow.ToString("0.000")));
                var fracHighs = String.Join(",", variantGroup.Select(mv => mv.FractionHigh.ToString("0.000")));
                var pValues = String.Join(",", variantGroup.Select(mv => ErrorProbabilityToQv(mv.PValue)));

                var vcfLine = String.Format(
                    "<{0}>\t{1}\t.\t{2}\t{3}\t.\tPASS\tCV={4};ZS={5};VF={6};CL={7};CH={8};PV={9}",
                    contig.Id, pos + 1, refBase, altBases, cov, counts, fracs, fracLows, fracHighs, pValues);

                if (contig.AminoVariants)
                {
                    var aaVariants = String.Join(",", variantGroup.Select(
                        mv => AminoVariant(contig.Sequence, pos, mv.AltBase)));

                    vcfLine = String.Format("{0};AA={1}", vcfLine, aaVariants);
                }

                _streamWriter.WriteLine(vcfLine);
            }
        }
    }
}
