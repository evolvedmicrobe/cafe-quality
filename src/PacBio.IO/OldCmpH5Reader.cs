using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.IO;
using System.Linq;
using PacBio.HDF;
using PacBio.Utils;
using PacBio.Align;


namespace PacBio.IO
{
    public interface IAlnSummary
    {
        uint AlnGroupID { get; }

        string MovieName { get; }

        uint RefGroup { get; }
        string ReferenceName { get; }

        int TemplateStart { get; }
        int TemplateEnd { get; }
        int TemplateLength { get; }

        Strand Strand { get; }
        int SubReadId { get; }
        int ReadStart { get; }
        int ReadEnd { get; }
        int ReadLength { get; }

        int Matches { get; }
        int Mismatches { get; }
        int Insertions { get; }
        int Deletions { get; }
        
        double Accuracy { get; }

        // Position of alignment in AlnArray
        int OffsetBegin { get; }
        int OffsetEnd { get; }

        string ReadName { get; }

        int HoleNumber { get; }

        string SequencingChemistry { get; }

        IAlignment ReadAlignment();
        Tuple<string, string> ReadAlignmentStrings(Orientation orientation = Orientation.Native);

        T[] ReadMetric<T>(string metricName);
    }


    public struct AlnSummary : IAlnSummary
    {
        private uint[] row;
        private OldCmpH5Reader cmpH5;
        public AlnSummary(uint[] row, OldCmpH5Reader cmp)
        {
            cmpH5 = cmp;
            this.row = row;
        }

        public uint AlnGroupID { get { return row[1]; } }

        public string MovieName { get { return cmpH5.MovieInfoIDMap[row[2]].Name; } }
        public string SequencingChemistry { get { return cmpH5.MovieInfoIDMap[row[2]].SequencingChemistry; } }

        public uint RefGroup { get { return row[3]; } }

        public string ReferenceName { get { return cmpH5.RefInfoIDMap[cmpH5.RefGroupIDMap[RefGroup].RefInfoID].FullName; } }

        public int TemplateStart { get { return (int)row[4]; } }
        public int TemplateEnd { get { return (int)row[5]; } }

        public int TemplateLength
        {
            get { return Math.Abs(TemplateEnd - TemplateStart); }
        }
        
        public Strand Strand
        {
            get { return row[6] == 0 ? Strand.Forward : Strand.Reverse; }
        }

        public int HoleNumber { get { return (int)row[7]; } }
        public int SetNumber { get { return (int)row[8]; } }
        public int StrobeNumber { get { return (int)row[9]; } }
        public int SubReadId { get { return (int)row[10]; } }

        public int ReadStart { get { return (int)row[11]; } }
        public int ReadEnd { get { return (int)row[12]; } }

        public int ReadLength
        {
            get { return ReadEnd - ReadStart;  }
        }

        // Map QV
        public int MapQV { get { return (int)row[13]; } }

        // Accuracy fields
        public int Matches { get { return (int)row[14]; } }
        public int Mismatches { get { return (int)row[15]; } }
        public int Insertions { get { return (int)row[16]; } }
        public int Deletions { get { return (int)row[17]; } }

        public double Accuracy
        {
            get { return 1.0 - (double)(Insertions + Mismatches + Deletions) / ReadLength; }
        }
        
        // Position of alignment in AlnArray
        public int OffsetBegin { get { return (int)row[18]; } }
        public int OffsetEnd { get { return (int)row[19]; } }

        public string ReadName { get { return String.Format("{0}/{1}", MovieName, HoleNumber); } }

        public IAlignment ReadAlignment()
        {
            return cmpH5.ReadAlignment(this);
        }

        public Tuple<string, string> ReadAlignmentStrings(Orientation orientation)
        {
            return cmpH5.ReadAlignmentStrings(this, orientation);
        }

        /// <summary>
        /// FIXME!! make a more type-safe / explicit list of the available metrics
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="metricName"></param>
        /// <returns></returns>
        public T[] ReadMetric<T>(string metricName)
        {
            return cmpH5.ReadMetric<T>(metricName, this);
        }
    }

    /// <summary>
    /// Routines for reading and translating HDF5 Alignment files (.cmp.h5).
    /// </summary>
    public class OldCmpH5Reader : IDisposable
    {

        //GLOBAL ALIGNMENT INDEX TABLE (SpringField)
        /*
        0 AlignmentId, 
        1 AlnGroupId, 
        2 MovieId,
        3 RefSeqId, 
        4 tStart, 
        5 tEnd, 
        6 RCRefStrand,
        7 HoleNumber, 
        8 SetNumber, 
        9 StrobeNumber, 
        10 SubreadId, 
        11 rStart, 
        12 rEnd,
        13 MapQV, 
        14 nM,
        15 nMM,
        16 nIns, 
        17 nDel, 
        18 offset_begin, 
        19 offset_end
        20 nBackRead  (optional, if sorted/indexed)
        21 nReadOverlap (optional, if sorted/indexed)
        */

        public class MovieInfoRecord
        {
            /// <summary>
            /// Movie ID in cmp.h5
            /// </summary>
            public uint ID { get; internal set; }

            /// <summary>
            /// Movie name
            /// </summary>
            public string Name { get; internal set; }

            /// <summary>
            /// Sequencing Chemistry
            /// </summary>
            public string SequencingChemistry { get; internal set; } 

            /// <summary>
            /// Camera frame rate
            /// </summary>
            public float FrameRate { get; internal set; }
        }

        public class RefInfoRecord
        {
            /// <summary>
            /// Reference ID
            /// </summary>
            public uint ID { get; internal set; }

            /// <summary>
            /// Full Name of contig
            /// </summary>
            public string FullName { get; internal set; }

            /// <summary>
            /// Size of contig
            /// </summary>
            public int Length { get; internal set; }

            /// <summary>
            /// MD5 hash of reference
            /// </summary>
            public string MD5 { get; internal set; }
        }
        
        public class RefGroupRecord
        {
            public uint ID { get; internal set; }
            public uint RefInfoID { get; internal set; }

            public string Path { get; internal set; }
        }

        public static OldCmpH5Reader CreateReader(string fileName)
        {
            IChunkFile fileHandle;
            try
            {
                fileHandle = HDFFile.Open(fileName, FileMode.Open, FileAccess.Read);
            }
            catch (Exception e)
            {
                string errMsg = e.Message;
                throw new IOException("Cannot open HDF file for reading: " + errMsg, e);
            }
            return new OldCmpH5Reader(fileHandle);
        }

        // Attribute if the cmp.h5 is sorted or not
        public bool IsSorted
        {
            get { return fileHandle.GetChild("/RefGroup/OffsetTable") != null; }
        }

        /// <summary>
        /// Direcotry that Primary Analysis results came from.
        /// </summary>
        public string ReportsFolder { get; private set; }


        // HDF5 file handle
        private IChunkFile fileHandle;

        // RefInfo
        public RefInfoRecord[] References;
        internal Dictionary<uint, RefInfoRecord> RefInfoIDMap;
        internal Dictionary<string, RefInfoRecord> RefInfoNameMap;
        internal Dictionary<string, RefInfoRecord> RefInfoMD5Map; 

        // MovieInfo
        public MovieInfoRecord[] Movies;
        internal Dictionary<uint, MovieInfoRecord> MovieInfoIDMap;
        internal Dictionary<string, MovieInfoRecord> MovieInfoNameMap;

        // RefGroup
        public RefGroupRecord[] RefGroups;
        internal Dictionary<uint, RefGroupRecord> RefGroupIDMap; 

        // AlnGroup
        internal Dictionary<uint, string> AlnGroupIdToPath;

        // Cache datasets to avoid metadata overhead
        private Dictionary<string, IDataset> datasetHandleCache = new Dictionary<string, IDataset>();

        private OldCmpH5Reader(IChunkFile fh)
        {
            var attr = fh.GetAttribute("ReportsFolder");
            if (attr != null)
            {
                ReportsFolder = attr.ReadSingleton<string>();
            }
            else
            {
                ReportsFolder = "Analysis_Results";
            }

            fileHandle = fh;
            Initialize();
        }


        private ReadOnlyCollection<IAlnSummary> alnCache = null;


        IDataset AlignmentIndex
        {
            get
            {
                return (IDataset)fileHandle.GetChild("AlnInfo/AlnIndex");
            }
        }


        /// <summary>
        /// An array of summary data structures of the alignments available.
        /// </summary>
        public IReadOnlyList<IAlnSummary> Alns
        {
            get
            {
                if (alnCache == null)
                {
                    var indexDs = AlignmentIndex;
                    var alnIndex = indexDs.Read() as uint[,];

                    var nHits = alnIndex.GetLength(0);
                    var nCols = alnIndex.GetLength(1);

                    var alns = new IAlnSummary[nHits];

                    for (int i = 0; i < nHits; i++)
                    {
                        var row = new uint[nCols];
                        for (int j = 0; j < nCols; j++)
                            row[j] = alnIndex[i, j];

                        alns[i] = new AlnSummary(row, this);
                    }

                    alnCache = new ReadOnlyCollection<IAlnSummary>(alns);
                }

                return alnCache;
            }
        }


        public IEnumerable<IAlnSummary> GetReadsForReference(uint refId)
        {
            if (!IsSorted)
                throw new ApplicationException(".cmp.h5 not sorted!");

            var offsets = ((IDataset)fileHandle.GetChild("/RefGroup/OffsetTable")).Read() as uint[,];

            var startEnd = Enumerable.Range(0, offsets.GetLength(0))
                .Where(i => offsets[i, 0] == refId)
                .Select(i => Tuple.Create(offsets[i, 1], offsets[i, 2]))
                .ToList();

            if (!startEnd.Any())
                return new IAlnSummary[0];

            var offStart = (int) startEnd.First().Item1;
            var offEnd = (int) startEnd.First().Item2;

            if (offEnd - offStart <= 0)
                return new IAlnSummary[0];

            return Enumerable.Range(offStart, offEnd - offStart)
                .Select(i => Alns[i]);
        }


        public IEnumerable<IAlnSummary> GetReadsInRange(uint refId, int refStart, int refEnd)
        {
            // TODO (lhepler) -- use bisection like David does in pbcore.io.CmpH5Reader
            return GetReadsForReference(refId)
                .Where(alnSummary => alnSummary.TemplateStart < refEnd &&
                                     alnSummary.TemplateEnd > refStart);
        }


        private Dictionary<string, Dictionary<int, List<IAlnSummary>>> movieHoleNumberCache;
        

        private void PrepMovieHoleNumberCache()
        {
            movieHoleNumberCache = new Dictionary<string, Dictionary<int, List<IAlnSummary>>>();

            foreach (var a in Alns)
            {
                Dictionary<int, List<IAlnSummary>> hnCache;
                if (!movieHoleNumberCache.TryGetValue(a.MovieName, out hnCache))
                {
                    hnCache = new Dictionary<int, List<IAlnSummary>>();
                    movieHoleNumberCache[a.MovieName] = hnCache;
                }

                List<IAlnSummary> alns;
                if (hnCache.TryGetValue(a.HoleNumber, out alns))
                {
                    alns.Add(a);
                }
                else
                {
                    hnCache[a.HoleNumber] = new List<IAlnSummary>() { a };
                }
            }
        }


        public List<IAlnSummary> GetAlignmentsForZmw(string movieName, int holeNumber)
        {
            if (movieHoleNumberCache == null)
            {
                PrepMovieHoleNumberCache();
            }

            var c1 = movieHoleNumberCache[movieName];
            List<IAlnSummary> res;
            if (c1.TryGetValue(holeNumber, out res))
            {
                return res;
            }
            else
            {
                return new List<IAlnSummary>();
            }
        }

        // ReSharper disable PossibleNullReferenceException
        private void MapReferences()
        {
            var refInfoID = GetDataset("RefInfo/ID").Read() as uint[];
            var refInfoName = GetDataset("RefInfo/FullName").Read() as string[];
            var refMD5 = GetDataset("RefInfo/MD5").Read() as string[];
            var refLength = GetDataset("RefInfo/Length").Read() as uint[];

            var nMovies = refInfoID.Length;

            References = nMovies.Fill(
                i =>
                new RefInfoRecord
                {
                    ID = refInfoID[i],
                    FullName = refInfoName[i],
                    Length = (int)refLength[i],
                    MD5 = refMD5[i]
                });

            RefInfoIDMap = References.ToDictionary(r => r.ID);
            RefInfoNameMap = References.ToDictionary(r => r.FullName);
            RefInfoMD5Map = References.ToDictionary(r => r.MD5);


            var refGroupID = GetDataset("RefGroup/ID").Read() as uint[];
            var refGroupRefInfoIDForeignKey = GetDataset("RefGroup/RefInfoID").Read() as uint[];
            var refGroupPath = GetDataset("RefGroup/Path").Read() as string[];

            RefGroups = refGroupID.Length.Fill(
                i => new RefGroupRecord
                {
                    ID = refGroupID[i],
                    RefInfoID = refGroupRefInfoIDForeignKey[i],
                    Path = refGroupPath[i]
                });

            RefGroupIDMap = RefGroups.ToDictionary(r => r.ID);
        }


        private void MapMovies()
        {
            // ID
            var dsId = GetDataset("MovieInfo/ID");
            var ids = dsId.Read() as uint[];
            var nMovies = ids.Length;

            // Name
            var dsName = GetDataset("MovieInfo/Name");
            var names = dsName.Read() as string[];

            // Frame rate w/ fallback
            var dsFrameRate = GetDataset("MovieInfo/FrameRate");
            var frameRate = nMovies.Fill(i => float.NaN);
            if (dsFrameRate != null)
                frameRate = dsFrameRate.Read() as float[];

            // Sequencing chemistry w/ fallback
            var unknown = nMovies.Fill<string>(i => null);

            var dsBindingKit = GetDataset("MovieInfo/BindingKit");
            var bindingKit = dsBindingKit != null ? dsBindingKit.Read() as string[] : unknown;

            var dsSequencingKit = GetDataset("MovieInfo/SequencingKit");
            var sequencingKit = dsSequencingKit != null ? dsSequencingKit.Read() as string[] : unknown;

            var dsSoftwareVersion = GetDataset("MovieInfo/SoftwareVersion");
            var softwareVersion = dsSoftwareVersion != null ? dsSoftwareVersion.Read() as string[] : unknown;

            var dsSequencingChemistry = GetDataset("MovieInfo/SequencingChemistry");
            var sequencingChemistry = dsSequencingChemistry != null ? dsSequencingChemistry.Read() as string[] : unknown;

            // Make movie
            Movies = nMovies.Fill(
                i =>
                {
                    string chem;

                    if (bindingKit[i] != null && sequencingKit[i] != null && softwareVersion[i] != null)
                        chem = PacBio.Data.Chemistry.DecodeTriple(bindingKit[i], sequencingKit[i], softwareVersion[i]);
                    else if (sequencingChemistry[i] != null)
                        chem = sequencingChemistry[i];
                    else
                        throw new PacBio.Data.ChemistryLookupException("Chemistry lookup information could not be found in cmp.H5!");

                    return new MovieInfoRecord
                    {
                        ID = ids[i],
                        Name = names[i],
                        FrameRate = frameRate[i],
                        SequencingChemistry = chem
                    };
                });

            MovieInfoIDMap = Movies.ToDictionary(m => m.ID);
            MovieInfoNameMap = Movies.ToDictionary(m => m.Name);
        }

        private void MapAlnGroups()
        {
            IDataset dsId = GetDataset("AlnGroup/ID");
            var ids = dsId.Read() as uint[];

            IDataset dsName = GetDataset("AlnGroup/Path");
            var names = dsName.Read() as string[];

            AlnGroupIdToPath = new Dictionary<uint, string>();

            var nGroups = ids.Length;
            for (int i = 0; i < nGroups; i++)
            {
                AlnGroupIdToPath[ids[i]] = names[i];
            }
        }
        // ReSharper restore PossibleNullReferenceException

        public string CommandLine
        {
            get
            {
                    var ds = (IDataset) fileHandle.GetChild("FileLog/CommandLine");
                    return ds.ReadSingleton<string>();
            }
        }


        public T[] ReadMetric<T>(string name, IAlnSummary aln)
        {
            string path = AlnGroupIdToPath[aln.AlnGroupID] + "/" + name;
            var ds = GetDataset(path);

            if (ds.Datatype.NativeType != typeof (T))
            {
                throw new ApplicationException(
                    String.Format("Incorrect datatype for metric: {0}.  Requested {1}, type in file is {2}", name, typeof(T), ds.Datatype.NativeType)
                    );
            }

            var size = aln.OffsetEnd - aln.OffsetBegin;
            var start = new long[] { aln.OffsetBegin, 0 };
            var stride = new long[] { 1, 1 };
            var count = new long[] { size, 1 };
            var block = new long[] { 1, 1 };
            
            var typedArr = new T[size];
            var arr = (Array) typedArr;

            var dSpace = ds.Dataspace;
            dSpace.SelectHyperslab(start, stride, count, block);
            ds.Read(ref arr, dSpace);

            return typedArr;
        }

        IDataset GetAlnArray(uint readGroupId)
        {
            string path = AlnGroupIdToPath[readGroupId] + "/" + "AlnArray";
            return GetDataset(path);
        }

        private long[] oneL = new long[]{ 1 } ;

        public Tuple<string, string> ReadAlignmentStrings(IAlnSummary aln, Orientation orientation = Orientation.Native)
        {
            // Read the alignment data from the AlnArray
            var alignmentDs = GetAlnArray(aln.AlnGroupID);
            var size = aln.OffsetEnd - aln.OffsetBegin;
            var start = new long[] { aln.OffsetBegin };   // Array with starting coordinates of hyperslab  
            var count = new long[] { size };              // Array specifying how many blocks to select from the dataspace, in each dimension 

            Array alignmentFrame = new byte[size];
            var alignmentDspace = alignmentDs.Dataspace;
            alignmentDspace.SelectHyperslab(start, oneL, count, oneL);
            alignmentDs.Read(ref alignmentFrame, alignmentDspace);

            // Unpack alignment data
            var doRc = aln.Strand == Strand.Reverse && orientation == Orientation.Genomic;
            var alignedStrings = PackedAlignment.GetAlignedStrings((byte[])alignmentFrame, doRc);

            // Handle any weirdness in cmp.h5 file
            if (alignedStrings.Item1.Length != alignedStrings.Item2.Length)
                throw new ApplicationException("Invalid cmp.h5");

            return alignedStrings;
        }

        public IAlignment ReadAlignment(IAlnSummary aln)
        {
            var alignStrings = ReadAlignmentStrings(aln, Orientation.Native);
            var alignedRef = alignStrings.Item1;
            var alignedRead = alignStrings.Item2;
            var cells = AlignCell.AlignCellsFromStrings(alignedRef, alignedRead, aln.TemplateStart, aln.ReadStart);
            var tplString = AlignCell.StripMarkup(alignedRef);
            var tpl = new TemplateSpec(tplString, aln.ReferenceName, aln.Strand, aln.TemplateStart);
            var readString = AlignCell.StripMarkup(alignedRead);

            return SimpleAlignment.OfAlignCellsAlignedRead(cells, tpl, readString);
        }

        public void Dispose()
        {
            Dispose(true);
            GC.SuppressFinalize(this);
        }

        ~OldCmpH5Reader()
        {
            Dispose(false);
        }

        protected virtual void Dispose(bool disposing)
        {
            if (disposing && fileHandle != null)
                fileHandle.Dispose();
        }

        private void Initialize()
        {
            MapMovies();
            MapReferences();
            MapAlnGroups();
        }


        private IDataset GetDataset(string path)
        {
            IDataset result;

            if (datasetHandleCache.TryGetValue(path, out result))
            {
                return result;
            }

            var ds = (IDataset)fileHandle.GetChild(path);

            if (ds != null)
                datasetHandleCache[path] = ds;

            return ds;
        }

    }

    /// <summary>
    /// Routines for storing and converting read-reference alignments
    /// </summary>
    public static class PackedAlignment
    {
        public static readonly char[] StrandInChar = { '+', '-' };
        
        private static char[][] byteToBases;
        private static Dictionary<char, char> baseComplement;

        static PackedAlignment()
        {
            byteToBases = CreateByteToBases(new [] { 'A', 'C', 'G', 'T', '-', 'N' });
            SetupComplement();
        }

        public static char Complement(char b)
        {
            switch (b)
            {
                case 'A': return 'T';
                case 'C': return 'G';
                case 'G': return 'C';
                case 'T': return 'A';
                case '-':
                case 'N':
                    return b;
                default:
                    throw new ArgumentException(
                        String.Format("cannot complement base {0}", b));
            }
        }

        public static Tuple<string, string> GetAlignedStrings(byte[] alignment, bool doRc)
        {
            var n = alignment.Length;
            var tpl = new char[alignment.Length];
            var read = new char[alignment.Length];

            for (int i = 0; i < n; i++)
            {
                var bp = byteToBases[alignment[doRc ? n - 1 - i : i]];
                tpl[i] = doRc ? Complement(bp[1]) : bp[1];
                read[i] = doRc ? Complement(bp[0]) : bp[0];
            }
            
            return new Tuple<string, string>(new string(tpl), new string(read));
        }


        private static char[][] CreateByteToBases(char[] alphabet)
        {
            var byteToBasesMap = new char[256][];
            for (int top = 0; top < 5; top++)
            {
                for (int bottom = 0; bottom <= 5; bottom++)
                {
                    //top = read ;bottom = reference
                    if (bottom == 5)
                    {
                        if (top < 4)
                            byteToBasesMap[(1 << (4 + top)) | 15] = new [] { alphabet[top], alphabet[bottom] };
                        else
                            byteToBasesMap[15] = new [] { alphabet[top], alphabet[bottom] };
                    }
                    else
                    {
                        if (top == 4 && bottom == 4)
                            byteToBasesMap[0] = new [] { alphabet[top], alphabet[bottom] };
                        else if (top == 4)
                            byteToBasesMap[((1 << (bottom)))] = new [] { alphabet[top], alphabet[bottom] };
                        else if (bottom == 4)
                            byteToBasesMap[(1 << top + 4)] = new []{ alphabet[top], alphabet[bottom] };
                        else
                            byteToBasesMap[((1 << (4 + top)) | 1 << bottom)] = new [] { alphabet[top], alphabet[bottom] };
                    }
                }
            }
            return byteToBasesMap;
        }

        private static void SetupComplement()
        {
            baseComplement = new Dictionary<char, char>();
            baseComplement['A'] = 'T';
            baseComplement['C'] = 'G';
            baseComplement['G'] = 'C';
            baseComplement['T'] = 'A';
            baseComplement['-'] = '-';
            baseComplement['N'] = 'N';
        }
    }
}
