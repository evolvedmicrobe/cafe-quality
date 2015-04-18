using System;
using System.Linq;
using System.IO;

namespace PacBio.IO.Fasta
{
    /// <summary>
    /// Simple class for writing FASTA-formatted output to a TextWriter or file.
    /// </summary>
    public class SimpleFASTQWriter : IFASTAWriter
    {
        private TextWriter writer;
        private int _lineWidth = 70;
        /// <summary>
        /// Fixed line width to use when formatting sequence
        /// </summary>
        public int LineWidth { get { return _lineWidth; } set { _lineWidth = value; } }

        // Stupid default constrcutor for MATLAB land
        public SimpleFASTQWriter()
        {
        }

        public SimpleFASTQWriter(string fileName)
        {
            writer = new StreamWriter(new BufferedStream(new FileStream(fileName, FileMode.Create)));
        }

        public SimpleFASTQWriter(TextWriter writer)
        {
            this.writer = writer;
        }

        public void WriteEntry(FASTQEntry entry)
        {
            writeEntry(writer, entry);
        }

        public void Open(string fileName)
        {
            if (writer != null)
                writer.Close();

            writer = new StreamWriter(new BufferedStream(new FileStream(fileName, FileMode.Create)));
        }

        public static void writeEntry(TextWriter writer, FASTQEntry entry)
        {
            writeEntry(writer, entry.Header, entry.GetSequence(), entry.GetQV(), 70);
        }

        public static void writeEntry(TextWriter writer, string name, string sequence, uint[] qv, int lineWidth)
        {
            lock (writer)
            {
                if (writer != null)
                {
                    writer.WriteLine("@" + name);
                    writer.WriteLine(sequence);
                    writer.WriteLine("+");

                    // Convert the integer quality
                    var qvstring =
                        new String(qv.Select(v => (Char) (Math.Min(93, v) + 33)).ToArray());

                    writer.WriteLine(qvstring);

                }
                else
                {
                    throw new Exception("FASTQWriter: you haven't opened a file to write to");
                }
            }
        }

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

        ~SimpleFASTQWriter()
        {
            Dispose(false);
        }

        public void WriteEntry(string name, string sequence)
        {
            throw new NotImplementedException();
        }
    }
}
