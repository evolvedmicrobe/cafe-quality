using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;

namespace PacBio.IO.Fasta
{
    /// <summary>
    /// Simple class for reading multi-FASTA files.  
    /// </summary>
    /// <remarks>
    ///   Reads one record at a time from the file,
    ///   keeping the file handle open.  Access to records is through the Enumerator.
    ///   <para>
    ///   Supports character-based filtering through the CharacterFilter property.
    ///   </para>
    ///   <para>
    ///   Not thread-safe.
    ///   </para>
    /// </remarks>
    public class SimpleFASTAReader : IDisposable, IEnumerable<FASTAEntry>
    {
        public static IEnumerable<FASTAEntry> ReadEntries(string fn)
        {
            using (var reader = new SimpleFASTAReader(fn))
            {
                foreach (var r in reader)
                    yield return r;
            }

        }

        /// <summary>
        /// Identifying character that starts the header of a FASTA entry
        /// </summary>
        private const char TAG_CHAR = '>';

        #region Enumerator definition

        private class SimpleFASTAEnumerator : IEnumerator<FASTAEntry>
        {
            private StringFASTAEntry current;
            private SimpleFASTAReader parent;

            public SimpleFASTAEnumerator(SimpleFASTAReader parent)
            {
                this.parent = parent;
            }

            #region IEnumerator<FASTAEntry> Members

            FASTAEntry IEnumerator<FASTAEntry>.Current
            {
                get { return current; }
            }

            #endregion

            #region IDisposable Members

            void IDisposable.Dispose()
            {
                // nothing to do here
            }

            #endregion

            #region IEnumerator Members

            object System.Collections.IEnumerator.Current
            {
                get { return current; }
            }

            bool System.Collections.IEnumerator.MoveNext()
            {
                current = parent.ReadNext();
                return current != null;
            }

            void System.Collections.IEnumerator.Reset()
            {
                parent.Reset();
            }

            #endregion
        }
        #endregion

        private string fileName;

        private bool forceLowerCase = false;
        public bool ForceLowerCase
        {
            get { return forceLowerCase; }
            set { forceLowerCase = value; }
        }

        private string currentHeaderLine = null;
		
		private StringBuilder builder;

        internal TextReader myReader = null;

        internal TextReader Reader
        {
            get
            {
                if (myReader == null)
                {
                    Stream myStream = new BufferedStream(new FileStream(fileName, FileMode.Open, FileAccess.Read));
                    myReader = new StreamReader(myStream);
                }
                return myReader;
            }
        }
		
        #region Character filtering

        public delegate bool SymbolFilter(char symbol);

        internal SymbolFilter characterFilter = null;

        public SymbolFilter CharacterFilter
        {
            get { return characterFilter; }
            set { characterFilter = value; }
        }

        internal string filterString(string sequence)
        {
            if (characterFilter == null) return sequence;
            char[] newChars = new char[sequence.Length];
            int i = 0;
            foreach (char c in sequence)
            {
                if (characterFilter(c))
                    newChars[i++] = c;
            }
            return new string(newChars, 0, i);
        }

        #endregion

        public SimpleFASTAReader(string fileName)
        {
            this.fileName = fileName;
        }


        public SimpleFASTAReader(TextReader reader)
        {
            this.myReader = reader;
        }

        /// <summary>
        /// Used by iterator
        /// </summary>
        internal void Reset()
        {
            if ( myReader!=null ) myReader.Close();
            myReader = null;
        }

        /// <summary>
        /// Callback used by enumerator
        /// </summary>
        /// <returns></returns>
        internal StringFASTAEntry ReadNext()
        {
            if (Reader.Peek()==-1)
                return null;
            string line;
            builder = new StringBuilder();
            String headerLine = currentHeaderLine;
            while ((line = Reader.ReadLine()) != null)
            {
                if (line.Length>0 && line[0] == TAG_CHAR)
                {
                    if (currentHeaderLine == null)
                    {
                        currentHeaderLine = (line.Length > 1 ? line.Substring(1) : "");
                        headerLine = currentHeaderLine;
                    }
                    else
                    {
                        currentHeaderLine = (line.Length > 1 ? line.Substring(1) : "");
                        break;
                    }
                }
                else
                {
                    builder.Append( filterString(line) );
                }
            }
            string sequence = builder.ToString();
            if (forceLowerCase)
                sequence = sequence.ToLower();
            StringFASTAEntry entry = new StringFASTAEntry(headerLine, sequence);
            return entry;
        }
		
		/// <summary>
        /// Read an initial chunk from reference file starting from position 0
        /// </summary>
		internal LazyFASTAEntry StreamReadNext(int refLen, int bufferThreshold)
        {	
            if (Reader.Peek()==-1)
                return null;
            string line;
			string fLine;
			int refEnd;
            builder = new StringBuilder();
            String headerLine = currentHeaderLine;
			
            while ((line = Reader.ReadLine()) != null)
            {
                if (line.Length>0 && line[0] == TAG_CHAR)
                {
                    if (currentHeaderLine == null)
                    {
                        currentHeaderLine = (line.Length > 1 ? line.Substring(1) : "");
                        headerLine = currentHeaderLine;
                    }
                    else
                    {
                        currentHeaderLine = (line.Length > 1 ? line.Substring(1) : "");
                        break;
                    }
                }
                else
                {
					fLine = filterString(line);
					builder.Append( fLine );
					if(builder.Length > bufferThreshold)
						break;
                }
			
            }
			refEnd = builder.Length - 1;
		    string sequence = builder.ToString();
            if (forceLowerCase)
                sequence = sequence.ToLower();
            LazyFASTAEntry entry = new LazyFASTAEntry(refLen, refEnd, headerLine, sequence, this);
            return entry;
        }
		
		/// <summary>
        /// Read more bases from reference file as necessary
        /// </summary>
		internal void StreamReadMore(int desiredEnd, LazyFASTAEntry entry)
		{
			string line;
			string fLine;
			
			while ((line = Reader.ReadLine()) != null)
			{
				fLine = filterString(line);
				builder.Append(fLine);
				if (desiredEnd - 1 <= builder.Length - 1)
					break;
			}
			entry.RefEnd = builder.Length - 1;
		    string sequence = builder.ToString();
            if (forceLowerCase)
                sequence = sequence.ToLower();
		    entry.Sequence = sequence;
		}
					
        /// <summary>
        /// Reads a single FASTA entry from this FASTA reader.  
        /// </summary>
        /// <remarks>
        /// Shouldn't be used in conjuction with iteration over the reader.
        /// </remarks>
        /// <returns>The entry or null if there are no entries</returns>
        public FASTAEntry ReadSingleEntry()
        {
            return ReadNext();
        }

		 public FASTAEntry StreamReadSingleEntry(int refLen)
        {
			int bufferThreshold = 50;
            return StreamReadNext(refLen, bufferThreshold);
        }
		
        public IEnumerator<FASTAEntry> GetEnumerator()
        {
            return new SimpleFASTAEnumerator(this);
        }
        
        System.Collections.IEnumerator System.Collections.IEnumerable.GetEnumerator()
        {
           throw new NotImplementedException(); 
        }

        #region Cleanup methods

        protected virtual void Dispose(bool disposing)
        {
            if (myReader != null)
            {
                myReader.Close();
                myReader = null;
            }
        }

        public void Dispose()
        {
            Dispose(true);
            GC.SuppressFinalize(this);
        }

        ~SimpleFASTAReader()
        {
            Dispose(false);
        }

        #endregion
    }
}
