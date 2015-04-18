using System;
using System.Runtime.InteropServices;

namespace PacBio.IO.Fasta
{
    public interface IFASTAWriter : IDisposable
    {
        void Open(string filename);
        void WriteEntry(string name, string sequence);
    }
}
