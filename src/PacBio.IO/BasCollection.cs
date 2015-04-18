using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using PacBio.Utils;

namespace PacBio.IO
{
    // "File of filenames"
    public static class Fofn 
    {
        public static string[] Filenames(string fofnFile)
        {
            var fofnDir = Path.GetDirectoryName(fofnFile);
            var fofnData = File.ReadAllLines(fofnFile)
                .Select(f => Path.IsPathRooted(f) ? f : Path.Combine(fofnDir, f))
                .ToArray();
            return fofnData;
        }
    }

    /// <summary>
    /// Provide access to reads from a collection of bas.h5 files
    /// </summary>
    public class BasCollection
    {
        private IBaseSource[] readers;
        private Dictionary<string, IBaseSource> movieReaders;
        private char[] baseMap;

        /// <summary>
        /// Open a BasCollection from a .fofn file
        /// </summary>
        /// <param name="fofnFile"></param>
        public static BasCollection FromFofn(string fofnFile)
        {
            var fofnData = Fofn.Filenames(fofnFile);
            try
            {
                return new BasCollection(fofnData);
            }
            catch (FileNotFoundException ex)
            {
                throw new FileNotFoundException(String.Format("Malformed fofn, member file not found: {0}", ex.FileName));
            }
        }

        /// <summary>
        /// Initialize a BasCollection from a list of bas.h5 files
        /// </summary>
        public BasCollection(string[] basFiles)
        {
            Init(basFiles);
        }

        private void Init(string[] basFiles)
        {
            readers = basFiles.Map(BaseReader.CreateSource);
            var groupReaders = readers.GroupBy(r => r.Movie.MovieName);

            movieReaders = new Dictionary<string, IBaseSource>();

            // The fofn may contain a combination of exploded .bax.h5 entries
            // and bas.h5 entries. Slurp up any multipart readers into a BaseMultiPartReader
            // and track as a single entry
            foreach (var group in groupReaders)
            {
                var myReaders = group.ToArray();
                if (myReaders.Length == 1)
                {
                    movieReaders[group.Key] = myReaders[0];
                }
                else
                {
                    movieReaders[group.Key] = new BaseMultiPartReader(myReaders.Map(v => v as BaseReader));
                }
            }

            baseMap = movieReaders.First().Value[0].Zmw.Movie.BaseMap;
        }

        /// <summary>
        /// Get the subreads for a given ZMW by ReadId
        /// </summary>
        public Subread[] GetSubreadsForZmw(string readId)
        {
            var parts = readId.Split('/');
            var movie = parts[0];
            var holeNumber = Int32.Parse(parts[1]);

            return GetSubreadsForZmw(movie, holeNumber);
        }

        /// <summary>
        /// Get the subreads for a given ZMW by movie, holeNumber
        /// </summary>
        /// <returns>The subreads for zmw.</returns>
        /// <param name="movie">Movie.</param>
        /// <param name="holeNumber">Hole number.</param>
        public Subread[] GetSubreadsForZmw(string movie, int holeNumber)
        {
            if (movieReaders.ContainsKey(movie))
            {
                var bases = movieReaders[movie].ByHoleNumber(holeNumber);
                return bases.Subreads();
            }
            else
            {
                var msg = String.Format("Movie ID was not found in input dataset: {0}. Please check inputs", movie);
                throw new ApplicationException(msg);
            }
        }

        /// <summary>
        /// Get the movie names of all the movie parts
        /// </summary>
        public string[] Movies
        {
            get { return movieReaders.Select(k => System.IO.Path.GetFileName(k.Key)).ToArray(); }
        }

        public char[] BaseMap
        {
            get { return baseMap; }
        }

        /// <summary>
        /// Get a subread aligned to the given IAlnSummary object
        /// </summary>
        public AlignedSubread GetSubread(IAlnSummary alignment)
        {
            var reader = movieReaders[alignment.MovieName];
            var bases = reader.ByHoleNumber(alignment.HoleNumber);

            var subreads = Subread.SubreadRegions(bases);

            var matchingSubread =
                subreads.FirstOrDefault(r => r.Start <= alignment.ReadStart && r.End >= alignment.ReadEnd);

            if (matchingSubread == null)
            {
                System.Diagnostics.Debugger.Break();
            }

            return new AlignedSubread(bases, matchingSubread, alignment);
        }
        
        /// <summary>
        /// Get a subread aligned to the given IAlnSummary object
        /// </summary>
        public IZmwBases GetRead(IAlnSummary alignment)
        {
            var reader = movieReaders[alignment.MovieName];
            var bases = reader.ByHoleNumber(alignment.HoleNumber);

            return bases;
        }

        /// <summary>
        /// Enumerate all HoleNumbers for a movie
        /// </summary>
        public IEnumerable<int> HoleNumber(string movie)
        {
            return movieReaders[movie].ZmwSource.ZmwIndexer.HoleNumber;
        }

        /// <summary>
        /// Enumerate all subreads
        /// </summary>
        public IEnumerable<Subread> Subreads
        {
            get { return readers.SelectMany(r => r).SelectMany(BasesHelpers.Subreads); }
        }

        /// <summary>
        /// Enumerate all reads
        /// </summary>
        public IEnumerable<IZmwBases> Reads
        {
            get { return readers.SelectMany(r => r); }
        }

        /// <summary>
        /// Total number of ZMWs in dataset
        /// </summary>
        public int NumZmws
        {
            get { return readers.Select(r => r.Count).Sum(); }
        }

        public IList<string> Chemistries
        {
            get { return readers.Select(reader => reader.Movie.SequencingChemistry).ToList(); }
        }
    }
}
