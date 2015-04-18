using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace PacBio.IO.Fasta
{
    /// <summary>
    /// Models a single entry from a FASTA file
    /// </summary>
    public abstract class FASTAEntry
    {
        protected string header;
        public string Header
        {
            get { return header; }
            set { header = value; }
        }

        internal FASTAEntry() { header = ""; }

        internal FASTAEntry(string header) { this.header = header; }

        public abstract string GetSequence();

        /// <summary>
        /// Returns sub-sequence from the sequence (forward) in this entry.
        /// </summary>
        /// <param name="start">0-based start of subsequence</param>
        /// <param name="end">Exclusive end of sequence</param>
        /// <returns></returns>
        public abstract string GetSubSequenceFwd(int start, int end);
		
		/// <summary>
        /// Returns sub-sequence from the reverse complement of the sequence in this entry.
        /// </summary>
        /// <param name="start">0-based start of subsequence</param>
        /// <param name="end">Exclusive end of sequence</param>
        /// <returns></returns>
		public string GetSubSequenceRev(int start, int end)
        {
			string seq = GetSubSequenceFwd(start,end);
			StringBuilder rcSeq =  new StringBuilder("");
			for (int i = seq.Length-1; i >=0 ; i--)
			{
				switch(seq[i])
				{
					case 'a':
					case 'A':
						rcSeq.Append('t');
					    break;
					
					case 'c':
					case 'C':
						rcSeq.Append('g');
					    break;
						
					case 't':
					case 'T':
						rcSeq.Append('a');
					    break;
						
					case 'g':
					case 'G':
						rcSeq.Append('c');
					    break;
					
					default:
						rcSeq.Append(seq[i]);
					    break;
				}				
			
			}
			return rcSeq.ToString();
        }
    }

    public class StringFASTAEntry : FASTAEntry
    {
        private string sequence;

        public StringFASTAEntry(string header, string sequence) : base( header )
        {
            this.sequence = sequence;
        }

        public override string GetSequence()
        {
            return sequence;
        }

        public override string GetSubSequenceFwd(int start, int end)
        {
            if (end < start)
                throw new NotSupportedException(string.Format(
                    "{0}<{1} is not allowed in GetSubSequence", end, start));
            if (end > sequence.Length)
                throw new IndexOutOfRangeException(string.Format("Trying to access index {0} from FASTA sequence (len={1})", 
                    end, sequence.Length) );
            return sequence.Substring(start, end - start);
        }
		
    }

    public class FASTQEntry : StringFASTAEntry
    {
        private uint[] qv;

        public FASTQEntry(string header, string sequence, uint[] qv) : base(header, sequence)
        {
            this.qv = qv;
        }

        public uint[] GetQV()
        {
            return qv;
        }
    }

    public class LazyFASTAEntry : FASTAEntry
    {   
		private SimpleFASTAReader reader;
		public int RefLen{get;set;}
		public int RefEnd{get;set;}
		public string Sequence{get;set;}

        public LazyFASTAEntry(int refLen, int refEnd,  string header, string sequence, SimpleFASTAReader reader) : base( header )
        {
            Sequence = sequence; //Portion of reference that is read into this entry
			RefEnd = refEnd; //Index of last character of reference read into this entry
			RefLen = refLen; //Length of entire reference
			this.reader = reader;
        }

		public override string GetSequence()
        {
            return Sequence;
        }
		
        public override string GetSubSequenceFwd(int start, int end)
        {
			if (end < start)
                throw new NotSupportedException(string.Format(
                    "{0}<{1} is not allowed in GetSubSequence", end, start));
            if (end > RefLen)
                throw new IndexOutOfRangeException(string.Format("Trying to access index {0} from FASTA sequence (len={1})", 
                    end, RefLen) );
			
			//Note here that end = exclusive and RefEnd = inclusive , in terms of indexing            
			if(end -1 > RefEnd)
			{
				reader.StreamReadMore(end, this);
			}
            return Sequence.Substring(start, end - start);
        }
		
    }
}
