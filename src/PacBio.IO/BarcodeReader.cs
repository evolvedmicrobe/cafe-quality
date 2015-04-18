using System.IO;
using System.Linq;
using System.Collections.Generic;
using PacBio.HDF;
using PacBio.Utils;

namespace PacBio.IO
{
    /// <summary>
    /// Read PacBio read barcoding files produced by the P_Barcode module of SMRT Pipe.
    /// Barcode results are typically stored in a series *.bc.h5 files, one file per movie-part, and a corresponding
    /// a barcode.fofn file
    /// </summary>
    public class BarcodeReader
    {
        /// <summary>
        /// Barcode call for a single ZMW
        /// </summary>
        public struct BarcodeCall
        {
            internal BarcodeCall(int row, MovieBarcodeReader reader)
            {
                this.row = row;
                this.reader = reader;
            }

            private readonly int row;
            private readonly MovieBarcodeReader reader;

            /// <summary>
            /// Movie containing this read
            /// </summary>
            public string Movie
            {
                get { return reader.Movie; }
            }

            /// <summary>
            /// Hole number of the read
            /// </summary>
            public int HoleNumber
            {
                get { return reader.Data[row, ColHoleNumber]; }
            }

            /// <summary>
            /// Read Id in MovieName/HoleNumber form
            /// </summary>
            public string ReadId
            {
                get { return Movie + @"/" + HoleNumber; }
            }

            /// <summary>
            /// Number of adapter hits observed in this read
            /// </summary>
            public int NumAdapters
            {
                get { return reader.Data[row, ColNumAdapters]; }
            }

            /// <summary>
            /// Barcode name of barcode called for this read
            /// </summary>
            public string BarcodeName
            {
                get { return reader.BarcodeNames[reader.Data[row, ColIndex1]]; }
            }

            /// <summary>
            /// Score of this barcode call
            /// </summary>
            public int Score
            {
                get { return reader.Data[row, ColScore1]; }
            }

            /// <summary>
            /// Barcode name of barcode called for this read
            /// </summary>
            public string BarcodeName2
            {
                get { return reader.BarcodeNames[reader.Data[row, ColIndex2]]; }
            }

            /// <summary>
            /// Score of this barcode call
            /// </summary>
            public int Score2
            {
                get { return reader.Data[row, ColScore2]; }
            }

            /// <summary>
            /// Average score of the barcode hits
            /// </summary>
            public int AvgScore
            {
                get { return (NumAdapters <= 0) ? 0 : Score / NumAdapters; }
            }

            /// <summary>>
            /// Barcode score ratio
            /// </summary>
            public double ScoreRatio
            {
                get { return (Score <= 0 || Score2 <= 0) ? 1.0 : ((double)Score) / Score2; }
            }
        }

        /// <summary>
        /// Initializes a multi-movie BarcodeReader from a barcode fofn file.
        /// </summary>
        public BarcodeReader(string fofnFile)
        {
            var fofnDir = Path.GetDirectoryName(fofnFile);
            Readers = File.ReadAllLines(fofnFile)
                .Select(f => Path.IsPathRooted(f) ? f : Path.Combine(fofnDir, f))
                .Select(f => new MovieBarcodeReader(f))
                .ToArray();
        }

        /// <summary>
        /// Set of movie part barcode readers
        /// </summary>
        public MovieBarcodeReader[] Readers
        {
            get;
            private set;
        }

        /// <summary>
        /// Enumerate all barcode calls made in this dataset
        /// </summary>
        public IEnumerable<BarcodeCall> Calls
        {
            get
            {
                return Readers.SelectMany(r => r.Calls);
            }
        }

        /// <summary>
        /// List of barccode names available in this dataset
        /// </summary>
        public string[] BarcodeNames
        {
            get
            {
                var bcSet = new HashSet<string>();
                Readers.ForEach(r => bcSet.UnionWith(r.BarcodeNames));
                return bcSet.ToArray();
            }
        }


        private const int ColHoleNumber = 0;
        private const int ColNumAdapters = 1;
        private const int ColIndex1 = 2;
        private const int ColScore1 = 3;
        private const int ColIndex2 = 4;
        private const int ColScore2 = 5;


        /// <summary>
        /// Class for reading barcoe calls from a single bc.h5 file
        /// </summary>
        public class MovieBarcodeReader
        {
            /// <summary>
            /// Construct a barcode reader for a single bc.h5 file
            /// </summary>
            /// <param name="barcodeFile"></param>
            public MovieBarcodeReader(string barcodeFile)
            {
                using (var f = HDFFile.Open(barcodeFile, System.IO.FileMode.Open, System.IO.FileAccess.Read))
                {
                    using (var best = (IDataset) ((IGroup) f.GetChild("BarcodeCalls")).GetChild("best"))
                    {

                        Movie = best.GetAttribute("movieName").ReadSingleton<string>();
                        BarcodeNames = (string[]) best.GetAttribute("barcodes").Read();

                        Data = (int[,]) best.Read();
                    }
                }
            }

            internal int[,] Data
            {
                get;
                set;
            }

            /// <summary>
            /// Barcode names available in this file
            /// </summary>
            public string[] BarcodeNames
            {
                get;
                private set;
            }

            /// <summary>
            /// Movie name of ZMWs in this barcode file
            /// </summary>
            public string Movie
            {
                get;
                private set;
            }

            /// <summary>
            /// Enumerate all barcode calls
            /// </summary>
            public IEnumerable<BarcodeCall> Calls
            {
                get
                {
                    for (int i = 0; i < Data.GetLength(0); i++)
                    {
                        yield return new BarcodeCall(i, this);
                    }
                }
            }
        }
    }
}

