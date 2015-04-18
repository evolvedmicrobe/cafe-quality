using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;

namespace PacBio.IO.Fasta
{
    /// <summary>
    /// Simple class for writing FASTA-formatted output to a TextWriter or file.
    /// </summary>
    public class SimpleFASTAWriter : IFASTAWriter
    {
        private TextWriter writer;
        private int _lineWidth = 70;
        /// <summary>
        /// Fixed line width to use when formatting sequence
        /// </summary>
        public int LineWidth { get { return _lineWidth; } set { _lineWidth = value; } }

        // Stupid default constrcutor for MATLAB land
        public SimpleFASTAWriter()
        {
        }

        public SimpleFASTAWriter(string fileName)
        {
            this.writer = new StreamWriter(new BufferedStream(new FileStream(fileName, FileMode.Create)));
        }

        public SimpleFASTAWriter(TextWriter writer)
        {
            this.writer = writer;
        }

        public void WriteEntry(string name, string sequence)
        {
            writeEntry( writer, name, sequence, LineWidth );
        }

        public void WriteEntry(FASTAEntry entry)
        {
            writeEntry(writer, entry);
        }

        public void Open(string fileName)
        {
            if (writer != null)
                writer.Close();

            this.writer = new StreamWriter(new BufferedStream(new FileStream(fileName, FileMode.Create)));
        }
		
		public static void writeEntry( TextWriter writer, FASTAEntry entry )
		{
			writeEntry( writer, entry.Header, entry.GetSequence(), 70 );
		}
		
		public static void writeEntry( TextWriter writer, string name, string sequence, int lineWidth )
		{
            lock (writer)
            {
                if (writer != null)
                {
                    writer.WriteLine(">" + name);
                    int len = sequence.Length;
                    for (int i = 0; i < len; i += lineWidth)
                    {
                        if (i + lineWidth <= len)
                            writer.WriteLine(sequence.Substring(i, lineWidth));
                        else
                            writer.WriteLine(sequence.Substring(i, len - i));
                    }
                }
                else
                {
                    throw new Exception("FastaWriter: you haven't opened a file to write to");
                }
            }
		}

        #region IDisposable Members

        protected virtual void Dispose(bool disposing)
        {
            if (writer != null)
            {
                writer.Close();
                writer = null;
            }
        }

        public void Dispose()
        {
            Dispose(true);
            GC.SuppressFinalize(this);
        }

        ~SimpleFASTAWriter()
        {
            Dispose(false);
        }

        #endregion
    }
}
