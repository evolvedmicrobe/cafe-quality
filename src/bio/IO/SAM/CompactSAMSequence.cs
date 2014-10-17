using System;
using System.Collections;
using System.Collections.Generic;
using System.Globalization;
using System.Text;
using System.Linq;
using Bio.Util;
using Bio.IO.SAM;

namespace Bio
{
    /// <summary>
    /// This class holds quality scores along with the sequence data.
    /// It is a container for the data in the SAMAlignedSequence, but much smaller to allow for faster parsing.
    /// much smaller.
    /// </summary>
    public class CompactSAMSequence : QualitativeSequence
    {
        #region SAMFields
        /// <summary>
        /// Name of Reference Sequence
        /// </summary>
        public string RName;
        /// <summary>
        /// Start position.
        /// </summary>
        public int Pos;

        /// <summary>
        /// SAM Flag vales
        /// </summary>
		public SAMFlags SAMFlags;

        /// <summary>
        /// The SAM Cigar Field
        /// </summary>
        public string CIGAR
		{
			get { return pCIGAR; }
			set {
				pCIGAR = value;
				if (!CigarUtils.NoInformationCigar (value)) {
					var alnLength = getRefSeqAlignmentLengthFromCIGAR ();
					this.RefEndPos = Pos + (alnLength > 0 ? alnLength - 1 : 0);
					if (RefEndPos < Pos) {
						throw new InvalidProgramException ();
					}
				}

			}
		}
		private string pCIGAR;
		/// <summary>
		/// Gets one based alignment end position of reference sequence depending on CIGAR Value.
		/// This value is INCLUSIVE!
		/// </summary>
		public int RefEndPos; 


        #endregion

		public CompactSAMSequence CreateTrimmedSequence(int newLength) {
			if (newLength > this.sequenceData.Length || newLength < 1)
			{
				throw new ArgumentOutOfRangeException("length");
			}

			byte[] newSequenceData = new byte[newLength];
			sbyte[] newQualityScores = new sbyte[newLength];

			Array.Copy (sequenceData, 0, newSequenceData, 0, newLength);
			Array.Copy (qualityScores, 0, newQualityScores, 0, newLength);
 
			//now to adjust cigar
			string newCigar = CIGAR;

			if (!CigarUtils.NoInformationCigar(CIGAR)) {
				var elements = CigarUtils.GetCigarElements (CIGAR);
				int curLen = 0;
				for (int i = 0; i < elements.Count; i++) {
					var ce = elements [i];
					var op = ce.Operation;
					if (op == CigarOperations.HARD_CLIP || op == CigarOperations.DELETION) {
						continue;
					} else if (op == CigarOperations.PADDING || op == CigarOperations.SKIPPED) {
						throw new NotImplementedException ();
					} else if (op == CigarOperations.SOFT_CLIP || op == CigarOperations.INSERTION || CigarOperations.OperationIsMatch (op)) {
						curLen += ce.Length;
						if (curLen == newLength) {
							newCigar = CigarUtils.CreateCigarString (elements.Take (i + 1));
							break;
						} else if (curLen > newLength) {
							var dif = curLen - newLength;
							ce.Length -= dif;
							elements [i]= ce;
							newCigar = CigarUtils.CreateCigarString (elements.Take (i + 1));
							break;
						}
					}
				}
			}
			var ns = new CompactSAMSequence (this.Alphabet, this.FormatType, newSequenceData, newQualityScores, false);
			ns.Pos = Pos;
			ns.CIGAR = newCigar;
			return ns;
		}

		public void TrimSequence(int newLength)
		{
			Count = newLength;
			if (newLength > this.sequenceData.Length || newLength <1)
			{
				throw new ArgumentOutOfRangeException("length");
			}
			Array.Resize (ref sequenceData, newLength);
			Array.Resize (ref qualityScores, newLength);

			//now to adjust cigar
			string newCigar = CIGAR;
			if (!CigarUtils.NoInformationCigar(CIGAR)) {
				var elements = CigarUtils.GetCigarElements (CIGAR);
				int curLen = 0;
				for (int i = 0; i < elements.Count; i++) {
					var ce = elements [i];
					var op = ce.Operation;
					if (op == CigarOperations.HARD_CLIP || op == CigarOperations.DELETION) {
						continue;
					} else if (op == CigarOperations.PADDING || op == CigarOperations.SKIPPED) {
						throw new NotImplementedException ();
					} else if (op == CigarOperations.SOFT_CLIP || op == CigarOperations.INSERTION || CigarOperations.OperationIsMatch (op)) {
						curLen += ce.Length;
						if (curLen == newLength) {
							newCigar = CigarUtils.CreateCigarString (elements.Take (i + 1));
							break;
						} else if (curLen > newLength) {
							ce.Length -= curLen - newLength;
							elements [i] = ce;
							newCigar = CigarUtils.CreateCigarString (elements.Take (i + 1));
							break;
						}
					}
				}
			}
			CIGAR = newCigar;
		}


		/// <summary>
		/// Gets the reference sequence alignment length depending on the CIGAR value.
		/// </summary>
		/// <returns>Length of the alignment.</returns>
		private int getRefSeqAlignmentLengthFromCIGAR()
		{
			if (CigarUtils.NoInformationCigar(CIGAR)) {
				return 0;
			}
			var elements = CigarUtils.GetCigarElements (CIGAR);
			int len = 0;
			foreach (var v in elements) {
				if (CigarUtils.CigarElementis_MDNX_Equal (v.Operation)) {
					len += v.Length;
				}
			}
			return len;
		}
        #region Constructors
        /// <summary>
        /// Initializes a new instance of the QualitativeSequence class with specified alphabet, quality score type,
        /// byte array representing symbols and encoded quality scores.
        /// Sequence and quality scores are validated with the specified alphabet and specified fastq format respectively.
        /// </summary>
        /// <param name="alphabet">Alphabet to which this instance should conform.</param>
        /// <param name="fastQFormatType">FastQ format type.</param>
        /// <param name="sequence">An array of bytes representing the symbols.</param>
        /// <param name="encodedQualityScores">An array of bytes representing the encoded quality scores.</param>
        public CompactSAMSequence(IAlphabet alphabet, FastQFormatType fastQFormatType, byte[] sequence, byte[] encodedQualityScores,bool validate)
            : base(alphabet, fastQFormatType, sequence, encodedQualityScores, validate)
        {
        }


        /// <summary>
        /// Initializes a new instance of the QualitativeSequence class with specified alphabet, quality score type,
        /// string representing symbols and encoded quality scores.
        /// Sequence and quality scores are validated with the specified alphabet and specified fastq format respectively.
        /// </summary>
        /// <param name="alphabet">Alphabet to which this instance should conform.</param>
        /// <param name="fastQFormatType">FastQ format type.</param>
        /// <param name="sequence">A string representing the symbols.</param>
        /// <param name="encodedQualityScores">A string representing the encoded quality scores.</param>
        public CompactSAMSequence(IAlphabet alphabet, FastQFormatType fastQFormatType, string sequence, string encodedQualityScores)
            : base(alphabet, fastQFormatType, sequence, encodedQualityScores, true)
        {
        }


        /// <summary>
        /// Initializes a new instance of the QualitativeSequence class with specified alphabet, quality score type,
        /// byte array representing symbols and signed byte array representing base quality scores 
        /// (Phred or Solexa base according to the FastQ format type).
        /// </summary>
        /// <param name="alphabet">Alphabet to which this instance should conform.</param>
        /// <param name="fastQFormatType">FastQ format type.</param>
        /// <param name="sequence">An array of bytes representing the symbols.</param>
        /// <param name="qualityScores">An array of signed bytes representing the base quality scores 
        /// (Phred or Solexa base according to the FastQ format type).</param>
        /// <param name="validate">If this flag is true then validation will be done to see whether the data is valid or not,
        /// else validation will be skipped.</param>
        public CompactSAMSequence(IAlphabet alphabet, FastQFormatType fastQFormatType, byte[] sequence, sbyte[] qualityScores, bool validate)
            : base(alphabet, fastQFormatType, sequence, qualityScores, validate) { }
        /// <summary>
        /// Initializes a new instance of the QualitativeSequence class with specified alphabet, quality score type,
        /// byte array representing symbols and integer array representing base quality scores 
        /// (Phred or Solexa base according to the FastQ format type).
        /// </summary>
        /// <param name="alphabet">Alphabet to which this instance should conform.</param>
        /// <param name="fastQFormatType">FastQ format type.</param>
        /// <param name="sequence">An array of bytes representing the symbols.</param>
        /// <param name="qualityScores">An array of integers representing the base quality scores 
        /// (Phred or Solexa base according to the FastQ format type).</param>
        /// <param name="validate">If this flag is true then validation will be done to see whether the data is valid or not,
        /// else validation will be skipped.</param>
        public CompactSAMSequence(IAlphabet alphabet, FastQFormatType fastQFormatType, byte[] sequence, int[] qualityScores, bool validate)
            : base(alphabet, fastQFormatType, sequence, qualityScores, validate) { }

        #endregion
	}
}
