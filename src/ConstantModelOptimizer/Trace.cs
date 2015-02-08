//using System;
//using System.Collections.Concurrent;
//using System.Collections.Generic;
//using System.IO;
//using System.Linq;
//using System.Threading;
//using System.Threading.Tasks;
//using PacBio.IO;
//using PacBio.Utils;
//using PacBio.Align;
//
//namespace ConstantModelOptimizer
//{
//    public class Trace
//    {
//        public TraceSet reader;
//        public string ReportsFolder { get{ return reader.ReportsFolder; } }
//        public string ReadID;
//        public IZmwBases ZmwBases { get; set; }
//        public IAlnSummary[] MultiAlignment { get; set; }
//        public IAlnSummary SmithWatermanAlignment { get; set; }
//        public IAlignment FullAlignment;
//
//        public float[] BaselineBias;
//        public float[] BaselineSigma;
//
//        public Trace(TraceSet reader, IZmwBases zmwBases, IAlnSummary alignment)
//        {
//            ReadID = ZmwBases.Zmw.Movie.MovieName + "/" + ZmwBases.Zmw.HoleNumber;
//            FullAlignment = reader.GetAlignment (alignment);
//            BaselineBias = reader.GetBaselineBias (ReadID);
//            BaselineSigma = reader.GetBaselineSigma (ReadID);
//
//        }
//
//        public int HoleNumber
//        {
//            get { return ZmwBases.Zmw.HoleNumber; }
//        }
//
//        public IMovieMetadata Metadata
//        {
//            get { return ZmwBases.Zmw.Movie;  }
//        }
//
//        public IAlignment GetAlignment(IAlnSummary s)
//        {
//            return Reader.ReadAlignment(s);
//        }
//    }
//
//}
//
