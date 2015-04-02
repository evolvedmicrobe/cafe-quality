using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using ConsensusCore;
using PacBio.Align;
using PacBio.IO;
using PacBio.IO.Fasta;
using PacBio.Utils;
using PacBio.HDF;
using System.ComponentModel;
using System.Reflection;

namespace PacBio.Consensus
{
    /// <summary>
    /// This enumeration lists all the ways that a ZMW could be processed by the CCS workflow.
    /// We return this value for each processed ZMW so that we can tally the different results an return them to the user.
    /// Each enumeration option is described by a description attribute which is used to produce the descrbipion in the text. 
    /// </summary>
    public enum CCSResultType { 
        [Description("Successful - Quiver consensus found")]
        Success = 0,
        [Description("Successful - But only 1 region, no true consensus")]
        OneRegionResult,
        [Description("Failed - Exception thrown")]
        Exception,
        [Description("Failed - ZMW was not productive")]
        NotProductive,
        [Description("Failed - Outside of SNR ranges")]
        OutsideSNR,
        [Description("Failed - No insert regions found")]
        NoInsertRegions,
        [Description("Failed - Not enough full passes")]
        NotEnoughFullPasses,
        [Description("Failed - Insert length too small")]
        InsertSizeTooSmall,
        [Description("Failed - Post POA requirements not met")]
        PostPOAFail,
        [Description("Failed - CCS Read below predicted accuracy")]
        PostCCSAccuracy,
        [Description("Failed - CCS Read was palindrome")]
        PostCCSPalindrome,
        [Description("Failed - CCS Read below SNR Threshold")]
        PostCCSSNR,
        [Description("Failed - CCS Read too short")]
        PostCCSShort
    }
    /// <summary>
    /// Keeps track of different CCS success and failures for all ZMWs.
    /// </summary>
    public class CcsResultsReport
    {
        /// <summary>
        /// An array to count the type of each result observed, elements will be filled as determined
        /// by the index in the CcsResult enum passed to the tally method.
        /// </summary>
        long[] resultsCount;

        /// <summary>
        /// Should a report be made? If not no information is stored or returned.
        /// </summary>
        bool makeReport;

        /// <summary>
        /// Initializes a new instance of the <see cref="PacBio.Consensus.CcsResultsReport"/> class.
        /// </summary>
        /// <param name="makeReport">If set to <c>true</c> make report.</param>
        public CcsResultsReport(bool makeReport)
        {
            this.makeReport = makeReport;
            if (makeReport)
            {
                int numOptions = Enum.GetNames(typeof(CCSResultType)).Length;
                resultsCount = new long[numOptions];
            }
        }

        /// <summary>
        /// Tallies the result in a thread safe way.
        /// </summary>
        /// <param name="result">Result.</param>
        public void TallyResult(CCSResultType result)
        {
            if (makeReport)
            {
                int index = (int)result;
                System.Threading.Interlocked.Increment(ref resultsCount[index]);
            }
        }

        /// <summary>
        /// The string we use to format the percentage of ZMWs in each group.
        /// </summary>
        internal readonly static string PERCENTAGE_FORMAT = "P2";

        /// <summary>
        /// Makes a final report based on all the results tallied.  If initialized 
        /// to not make a report, simply returns an empty string.
        /// </summary>
        /// <returns>The final report.</returns>
        public string MakeFinalReport()
        {
            if (!makeReport)
            {
                return String.Empty;
            }
            double totalZMWs = (double)resultsCount.Sum();
            // Make a line for each result
            var results = new List<string[]>(resultsCount.Length + 1);
            results.Add( new string[] {"Zmw Result" , "#-Zmws", "%-Zmws"});
            for (int i = 0; i < resultsCount.Length; i++)
            {
                var description = getCcsResultDescription(i);
                var cnt = resultsCount[i];//.ToString();
                var perc = (cnt / totalZMWs).ToString(PERCENTAGE_FORMAT);
                results.Add(new string[] { description, cnt.ToString(), perc });
            }

            // Format each line by adding white space
            int bufferCharsPerCol = 5; // how many extra characters to add.
            for (int i = 0; i < results.First().Length; i++)
            {
                var maxSize = results.Max(z => z[i].Length);
                int newSize = maxSize + bufferCharsPerCol;
                foreach(var p in results) {
                    p[i] = p[i] + new string(' ', newSize - p[i].Length);
                }
            }

            // Now combine into one big string
            System.Text.StringBuilder sb = new System.Text.StringBuilder();
            sb.AppendLine("Result Report for the " + ((long)totalZMWs) + " Zmws processed");  
            foreach (var res in results)
            {
                sb.AppendLine(String.Join("", res));
            }
            return sb.ToString();
        }


        /// <summary>
        /// Get the description for a CcsResult enum given a value.
        /// </summary>
        /// <returns>The ccs result description.</returns>
        /// <param name="value">The integer cast version of the enum.</param>
        private static string getCcsResultDescription(int value)
        {
            var converted = (CCSResultType)value;
            FieldInfo fi = typeof(CCSResultType).GetField(converted.ToString());
            DescriptionAttribute[] attributes =
                (DescriptionAttribute[])fi.GetCustomAttributes(
                    typeof(DescriptionAttribute),
                    false);

            if (attributes != null &&
                attributes.Length > 0)
                return attributes[0].Description;
            else
                throw new NotImplementedException("Could not acquire description of CCS result type. "
                   + "This attribute likely needs to be added to the CcsResult enum.");
        }

    }



    public class CsvSink : BasicSinkStage<Tuple<IZmwBases, IZmwConsensusBases>>
    {
        private CsvWriter csvWriter;

        public CsvSink(string filename)
        {
            if (!String.IsNullOrEmpty(filename))
            {
                csvWriter = new CsvWriter(filename);
            }
        }

        public override void OnNext(Tuple<IZmwBases, IZmwConsensusBases> reads)
        {
            var rawRead = reads.Item1;
            var ccsRead = reads.Item2;

            if (csvWriter != null  && 
                ccsRead.Zmw.ZmwType == ZmwType.Sequencing &&
                ccsRead.Sequence.Length > 0)
            {
                csvWriter.AddCol("MovieName", ccsRead.Zmw.Movie.MovieName);
                csvWriter.AddCol("HoleNumber", ccsRead.Zmw.HoleNumber);
                csvWriter.AddCol("PredictedRawAccuracy", rawRead.Metrics.ReadScore);
                csvWriter.AddCol("SnrT", rawRead.Metrics.HQRegionSNR[0]);
                csvWriter.AddCol("SnrG", rawRead.Metrics.HQRegionSNR[1]);
                csvWriter.AddCol("SnrA", rawRead.Metrics.HQRegionSNR[2]);
                csvWriter.AddCol("SnrC", rawRead.Metrics.HQRegionSNR[3]);
                csvWriter.AddCol("CCSReadLength", ccsRead.NumBases);
                csvWriter.AddCol("NumPasses", ccsRead.NumPasses);
                csvWriter.AddCol("PredictedCCSAccuracy", ccsRead.PredictedAccuracy);
                csvWriter.Row();
            }
        }

        public override void OnComplete()
        {
            if (csvWriter != null)
                csvWriter.Dispose();
        }

        public override void OnError(Exception e)
        {
            if (csvWriter != null)
                csvWriter.Dispose();
        }
    }


    /// <summary>
    /// Class for streaming FASTA / FASTQ writing during pipeline execution
    /// </summary>
    /// <typeparam name="T"></typeparam>
    public class FastaSink<T> : BasicSinkStage<T> where T : IZmwBases
    {
        /// <summary>
        /// Setup with MovieMetadata so that we know what ReadId convention to use.
        /// Set the fastaFile or fastqFile uris to null to prevent that file from being generated.
        /// </summary>
        public FastaSink(IMovieMetadata metadata, string fastaFile, string fastqFile)
        {
            if (!String.IsNullOrEmpty(fastaFile))
            {
                Log(LogLevel.INFO, "Opening FASTA output file '{0}'", fastaFile);
                var file = TempFile(fastaFile);
                fastaWriter = new SimpleFASTAWriter(file);
            }

            if (!String.IsNullOrEmpty(fastqFile))
            {
                Log(LogLevel.INFO, "Opening FASTQ output file '{0}'", fastqFile);
                var file = TempFile(fastqFile);
                fastqWriter = new SimpleFASTQWriter(file);
            }

            fastaTag = TraceReference.CreateSpringfield;

            // If we are a CCS read, then use the ccs read name
            if (typeof (IZmwConsensusBases).IsAssignableFrom(typeof (T)))
            {
                fastaTag = TraceReference.CreateSpringfieldCCS;
            }
            
            subreadFastaTag = TraceReference.CreateSpringfieldSubread;
        }

        private readonly Func<IZmwBases, string> fastaTag;
        private bool writeSequencingZmwsOnly = true;

        // tag = Func(IZmwBases, start, end)
        private Func<IZmwBases, int, int, string> subreadFastaTag;

        private SimpleFASTAWriter fastaWriter;
        private SimpleFASTAWriter subreadFastaWriter;
        private SimpleFASTQWriter fastqWriter;
        private SimpleFASTQWriter subreadFastqWriter;

        /// <summary>
        /// Provide this to re-apply a 'Productivity' length cut on adapter-trimmed subreads. 
        /// </summary>
        public int MinSubreadLength
        {
            get { return 35; }
        }

        // Only write Sequencing Zmws
        bool DoWrite(IZmwBases b)
        {
            if (writeSequencingZmwsOnly)
            {
                return b.Zmw.ZmwType == ZmwType.Sequencing;
            }

            return true;
        }


        /// <summary>
        /// Write an IZmwBases entry to the FASTA file
        /// </summary>
        private void WriteToFasta(IZmwBases read)
        {
            var sequence = read.Sequence;
            if (sequence.Length > 0)
            {
                var tag = fastaTag(read);
                fastaWriter.WriteEntry(tag, sequence);
            }
        }

        /// <summary>
        /// Write an IZmwBases entry to the FASTQ file
        /// </summary>
        /// <param name="read"></param>
        private void WriteToFastq(IZmwBases read)
        {
            var sequence = read.Sequence;
            if (sequence.Length > 0)
            {
                var qv = read.QV.Select(v => (uint)v).ToArray();
                var tag = fastaTag(read);

                var entry = new FASTQEntry(tag, sequence, qv);

                fastqWriter.WriteEntry(entry);
            }
        }

        /// <summary>
        /// Write subreads from an IZmwBases entry to the subreads FASTA and/or FASTQ file(s).
        /// </summary>
        /// <param name="read"></param>
        private void WriteToSubreadFastaq(IZmwBases read)
        {
            var productivity = read.Metrics.Productivity;
            var sequence = read.Sequence;
            var minLength = MinSubreadLength;

            if (productivity == ProductivityClass.Productive && sequence.Length > 0)
            {
                // Should never be null, since Prod==1
                var hqr = read.HQRegion();
                var hqrStart = hqr.Start;
                var hqrEnd = hqr.End;

                var qvs = read.QV.Select(v => (uint)v).ToArray();

                // Get all the subreads
                foreach (var reg in read.Metrics.Regions)
                {
                    if (!reg.Type.Equals(RegionAnnotator.InsertRegionType))
                        continue;

                    // Insert region(s) may be outside the HQ region:
                    // Defined as the complement of the adapter regions over the read.
                    //
                    var start = Math.Max(reg.Start, hqrStart);
                    var end = Math.Min(reg.End, hqrEnd);
                    var len = end - start;

                    // Skip any really short stuff
                    if (len < minLength)
                        continue;

                    var seq = sequence.Substring(start, len);
                    var tag = subreadFastaTag(read, start, end);

                    if (subreadFastaWriter != null)
                    {
                        subreadFastaWriter.WriteEntry(tag, seq);
                    }

                    if (subreadFastqWriter != null)
                    {
                        var qv = qvs.Skip(start).Take(len).ToArray();
                        var entry = new FASTQEntry(tag, seq, qv);
                        subreadFastqWriter.WriteEntry(entry);
                    }
                }
            }
        }

        public override void OnNext(T read)
        {
            // Skip this read if it doesn't pass the filter
            if (!DoWrite(read))
                return;

            // Send to the writers if they were constructed

            if (fastaWriter != null)
                WriteToFasta(read);

            if (fastqWriter != null)
                WriteToFastq(read);

            if (subreadFastaWriter != null || subreadFastqWriter != null)
                WriteToSubreadFastaq(read);
        }


        private void CloseWriters()
        {
            // Close the fasta files.

            if (fastaWriter != null)
            {
                fastaWriter.Dispose();
                fastaWriter = null;
            }

            if (subreadFastaWriter != null)
            {
                subreadFastaWriter.Dispose();
                subreadFastaWriter = null;
            }

            if (fastqWriter != null)
            {
                fastqWriter.Dispose();
                fastqWriter = null;
            }

            if (subreadFastqWriter != null)
            {
                subreadFastqWriter.Dispose();
                subreadFastqWriter = null;
            }
        }

        public override void OnComplete()
        {
            CloseWriters();
            FinalizeTempFiles();
        }

        public override void OnError(Exception e)
        {
            CloseWriters();
            DeleteTempFiles();
        }
    }

    public class ConsensusBaseSink : BasicSinkStage<IZmwConsensusBases>
    {
        private IChunkFile baseFile;
        private ConsensusBaseWriter baseWriter;

        /// <summary>
        /// Setup a ConsensusBaseWriter
        /// </summary>
        /// <param name="baseFileName">Target HDF5 file</param>
        /// <param name="scanDataGroup">Name of group to write to</param>
        /// <param name="hdfGroupName"></param>
        public ConsensusBaseSink(string baseFileName, IGroup scanDataGroup, string hdfGroupName)
        {
            Log(LogLevel.INFO, "Opening HDF5 output file '{0}'", baseFileName);

            baseFile = HDFFile.Open(TempFile(baseFileName), FileMode.OpenOrCreate, FileAccess.ReadWrite);
            baseWriter = new ConsensusBaseWriter(baseFile, hdfGroupName);

            // Copy ScanData from Pulse file to Base file
            var scanData = scanDataGroup as HDFGroup;
            if (scanData == null)
            {
                var msg = "An invalid ScanData group was passed to BaseSink";
                Log(LogLevel.ERROR, msg);
                throw new Exception(msg);
            }

            // Make sure that there isn't already a ScanData group there.
            baseFile.RemoveChild("ScanData");
            scanData.Copy(baseFile);
        }

        public void AddChemistryInformation(string sequencingChemistry)
        {
            baseWriter.AddChemistryInformation(sequencingChemistry);
        }

        public void AddChemistryInformation(string bindingKit, string sequencingKit, string changelistId)
        {
            baseWriter.AddChemistryInformation(bindingKit, sequencingKit, changelistId);
        }

        private void CloseFiles()
        {
            baseWriter.Dispose();
            baseWriter = null;

            baseFile.Dispose();
            baseFile = null;
        }

        public override void OnNext(IZmwConsensusBases item)
        {
            baseWriter.BufferedWrite(item, 100);
        }

        public override void OnComplete()
        {
            CloseFiles();
            FinalizeTempFiles();
        }

        public override void OnError(Exception e)
        {
            CloseFiles();
            DeleteTempFiles();
        }
    }

    public class ConsensusConfig
    {

        /// <summary>
        /// Minimum number of full passes for CCS reads. Must replace default value
        /// </summary>
        public int MinFullPasses = int.MaxValue;

        /// <summary>
        /// Minimum predicted accuracy for CCS reads. Must replace default value.
        /// </summary>
        public float MinPredictedAccuracy = 1.0f;

        /// <summary>
        /// Minimum length of CCS reads.
        /// </summary>
        public int MinLength = 1;

        /// <summary>
        /// Maximum length of CCS reads.
        /// </summary>
        public int MaxLength = int.MaxValue;

        public SnrCut SnrCut = SnrCut.PassAll;
    }


    /// <summary>
    /// Pipeline stage implementing CircularConsensus Sequencing.  For each ZMW the inputs are IZmwPulses, and IZmwBases.
    /// The output is IZmwConsensusBases
    /// </summary>
    public class CCSStream : PipelineMapper<ConsensusConfig, IZmwBases, Tuple<CCSResultType, IZmwConsensusBases>>
    {
        /// <summary>
        /// Instantiate the CCS algo with standard parameters
        /// </summary>
        public static CCSStream DefaultConfig
        {
            get
            {
                var cfg = new ConsensusConfig
                {
                    MinFullPasses = 2,
                    MinPredictedAccuracy = 0.9f
                };

                return new CCSStream(cfg);
            }
        }

        /// <summary>
        /// Instantiate the CCS algo with standard parameters
        /// </summary>
        public static CCSStream TrainingConfig
        {
            get
            {
                var cfg = new ConsensusConfig
                {
                    MinFullPasses = 0,
                    MinPredictedAccuracy = 0.1f
                };

                return new CCSStream(cfg);
            }
        }


        public int AdapterPadBases { get; private set; }
        public string Adapter { get; private set; }

        private ReadPartition partitioner;



        public CCSStream(ConsensusConfig config)
        {
            Config = config;

            Log(LogLevel.INFO, "Initializing Circular Consensus...");
            Log(LogLevel.INFO, "CCS Filter Settings -- MinPredictedAccuracy: {0}, MinFullPasses: {1}",
                Config.MinPredictedAccuracy, Config.MinFullPasses);
            //var algo = RecursionAlgo.Viterbi;

            // FIXME -- pipe in adapter
            Adapter = "ATCTCTCTCttttcctcctcctccgttgttgttgttGAGAGAGAT".ReverseComplement();
            partitioner = new ReadPartition(Adapter);
            AdapterPadBases = 8;

        }



        public override Tuple<CCSResultType, IZmwConsensusBases> Map(IZmwBases bases)
        {
            var zmw = bases.Zmw;
            Log(LogLevel.DEBUG, "CCS processing '{0}/{1}'", zmw.Movie.MovieName, zmw.HoleNumber);
            Tuple<CCSResultType, IZmwConsensusBases> toReturn;
            try
            {
                toReturn = InnerMap(bases);
            }
            catch (Exception e)
            {
                Log(LogLevel.ERROR, "CCS failed processing: '{0}/{1}/{2}'.  Hole Number: '{3}'", zmw.Movie.MovieName,
                    zmw.X,
                    zmw.Y, zmw.HoleNumber);
                Log(LogLevel.ERROR, "{0} caught during trace processing: {1}. Returning empty result", e.GetType(),
                    e.Message);
                Log(LogLevel.ERROR, "Stack Trace: {0}", e.StackTrace);

                toReturn = new Tuple<CCSResultType, IZmwConsensusBases>(CCSResultType.Exception, ZmwConsensusBases.Null(zmw));
            }
            Log(LogLevel.DEBUG, "CCS Result: {0},  '{1}/{2}'", toReturn.Item1 , zmw.Movie.MovieName, zmw.HoleNumber);
            return toReturn;
        }

        /// <summary>
        /// Return the longest region with a predicted accuracy > 0.75
        /// </summary>
        public ZmwConsensusBases PassBestRegion(IZmwBases bases, IList<DelimitedSeqReg> regions)
        {
            Func<DelimitedSeqReg, float> getAccPred = r =>
            {

                // Get the predicted accuray of the return region only.
                var qvs = bases.QV.Slice(r.Start, r.Length);
                return 1.0f - (float) qvs.Select(QVs.PhredProb).Average();
            };

            var bestRegion = regions.OrderByDescending(r => r.Length).FirstOrDefault(r => getAccPred(r) > Config.MinPredictedAccuracy);

            if (bestRegion == null)
            {
                return ZmwConsensusBases.Null(bases.Zmw);
            }
            else
            {
                var predictedAccuracy = getAccPred(bestRegion);
                return ZmwConsensusBases.PassRegion(bases.Zmw, bases, bestRegion, predictedAccuracy, bestRegion.Length);
            }
        }


       float RecalibratedAccuracyPrediction(string sequencingChemistry, float naivePrediction, int numPasses)
        {
            // Analysis was done here:
            //   /home/UNIXHOME/dalexander/Projects/Analysis/CCS-Performance/Calibration3   (P4-C2)
            //   /home/UNIXHOME/dalexander/Projects/Trainings/P6-C4/CCS/P5-C3-Recalibration (P5-C3)
            //   /home/UNIXHOME/dalexander/Projects/Trainings/P6-C4/CCS/DebuggingP6/Revisit-MH2155 (P6-C4)
            //
            // I essentially just fitted a linear model and then
            // subtracted one from the offset to guarantee a bit of
            // conservatism...
            //
            // This approach is all kinds of unsustainable and
            // horrible---when we get a better understanding of the
            // model we should fix the calibration problem at the
            // root.
            float predictedQSlope, numPassesSlope, interactionSlope, offset;

            // (Note: for one-pass "CCS" goes through PassBestRegion and seems reasonably calibrated; otherwise we use this
            //  linear model recalibration)
            if (sequencingChemistry == "P4-C2" && numPasses > 1)
            {
                predictedQSlope = 0.8569f;
                numPassesSlope = -1.6318f;
                interactionSlope = 0.0407f;
                offset = 7.7446f;
            } else if (sequencingChemistry == "P5-C3" && numPasses > 1) {
                predictedQSlope = 0.9022f;
                numPassesSlope = -0.7328f;
                interactionSlope = 0.0102f;
                offset = 6.0696f;
            } else if (sequencingChemistry == "P6-C4" && numPasses > 1) {
                predictedQSlope = 0.93610f;
                numPassesSlope = -0.97291f;
                interactionSlope = 0.02105f;
                offset = 5.54577f;
            } else {
                return naivePrediction;
            }

            float predictedQ = QVs.ProbToQV(1 - naivePrediction);
            float rawRecalibratedQ =
                (predictedQ * predictedQSlope +
                 numPasses * numPassesSlope +
                 predictedQ * numPasses * interactionSlope +
                 offset - 1);
            double recalibratedQ = Math.Min(30, Math.Max(0, rawRecalibratedQ));  // Clip to [0, 30]
            return (float)(1.0 - QVs.PhredProb(recalibratedQ));
        }

        /// <summary>
        /// Compute SMC for one trace.
        /// </summary>
        public Tuple<CCSResultType, IZmwConsensusBases> InnerMap(IZmwBases bases)
        {
            // Don't even bother if this is not a productive ZMW. 
            // This will filter out low SNR or otherwise crappy, which take a lot of time
            if (bases.Metrics.Productivity != ProductivityClass.Productive) // || bases.Zmw.HoleNumber != 80745)
                return new Tuple<CCSResultType, IZmwConsensusBases> (CCSResultType.NotProductive, ZmwConsensusBases.Null(bases.Zmw));


            if (!Config.SnrCut.Contains(bases.Metrics.HQRegionSNR))
                return new Tuple<CCSResultType, IZmwConsensusBases>(CCSResultType.OutsideSNR, ZmwConsensusBases.Null(bases.Zmw));

            //perf.Time(zmw.HoleNumber, "CCSPartition");

            // Find adapters
            var readRegions = partitioner.GetPartition(bases, AdapterPadBases);

            // 0 inserts -- empty result
            if (readRegions.InsertRegions.Count == 0)
                return new Tuple<CCSResultType, IZmwConsensusBases>(CCSResultType.NoInsertRegions, ZmwConsensusBases.Null(bases.Zmw));

            // In 2.0 we only emit >0-length sequences for true 'CCS' reads -- best estimate reads & CCS are moving to secondary
            // However, we supply the metrics for 'Read of Insert' reporting here, via the metrics attached to IZmwConsensusBases.
            // If there is not 1 full pass, just report the stats, not the sequence of the longest subread.
            // We will run CCS on everything else.  If it dones't make the 90% threshold, then we will not report any sequence for it, but we will report the metics.

            if (readRegions.InsertRegions.Count(r => r.AdapterHitAfter && r.AdapterHitBefore) < Config.MinFullPasses)
            {
                return new Tuple<CCSResultType, IZmwConsensusBases>(CCSResultType.NotEnoughFullPasses, ZmwConsensusBases.Null(bases.Zmw));
            }

            // Compute insert length
            var insertLength = readRegions.InsertRegions.Average(r => r.Length);
            //perf.Measure(zmw.HoleNumber, "AverageInsert", insertLength);

            // If this is an adapter dimer, or a very short insert then bail
            if (insertLength < 10)
                return new Tuple<CCSResultType, IZmwConsensusBases>(CCSResultType.InsertSizeTooSmall, ZmwConsensusBases.Null(bases.Zmw));


            //perf.Time(zmw.HoleNumber, "ComputePOA");

            // Quality metrics of the POA consensus
            float graphScore = 0.0f;
            //List<MutationScore> graphMutations = null;

            // The POA step must emit an initial template, and the subread regions to use
            TrialTemplate initialTpl = null;
            List<AlignedSequenceReg> regions = null;
            //List<MutationScore> startMutations = null;

            //var numCompletePasses = readRegions.InsertRegions.Count(r => r.AdapterHitBefore && r.AdapterHitAfter);

            // New-style POA setup -- use in all cases -- do a local alignment between subreads to find
            // overlapping sections
            initialTpl = FindConsensus.InitialConsensusTemplate(readRegions, bases, out graphScore,
                                                                out regions, AdapterPadBases, Adapter);

            // let the actual sequence vary by up to 5% of its estimable length
            // FIXME (lhepler) -- figure out the actual distribution
            // of differences in POA vs Quiver length, use 99th percentile here
            int tplVar = (int)Math.Ceiling(0.05 * initialTpl.Sequence.Length);
            int minTplLength = Math.Max(1, Config.MinLength - tplVar);
            // don't overflow
            int maxTplLength = (int.MaxValue - tplVar < Config.MaxLength) ? int.MaxValue : Config.MaxLength + tplVar;

            // Check that we have enough passes in play after the POA has finished
            if (regions.Count < Config.MinFullPasses ||
                initialTpl.Sequence.Length < minTplLength ||
                initialTpl.Sequence.Length > maxTplLength)
            {
                return new Tuple<CCSResultType, IZmwConsensusBases>(CCSResultType.PostPOAFail,  ZmwConsensusBases.Null(bases.Zmw));
            }

            IZmwConsensusBases result;

            // Local POA may return no regions -- just return the longest original region
            if (regions.Count < 2)
            {
                return new Tuple<CCSResultType, IZmwConsensusBases>(CCSResultType.OneRegionResult, PassBestRegion(bases, readRegions.InsertRegions));
            }
            else
            {
                // Refine the initial guess using a more sophisticated scoring scheme
                int maxIterations = 7;
                int iterationsTaken = 0; // in order T, G, A, C
                var snrs = bases.Metrics.HQRegionSNR;
                using (var snr = new SNR (snrs [2], snrs [3], snrs [1], snrs [0])) {
                    using (var ctx_params = new ContextParameters (snr)) {
                       /* Lance had set these to hard code 18, -12.5 if chemistry known
                        * Otherwise 24, -50 if chemistry is unspecified.  Not sure which is 
                        * correct, going to use 18, -12.5 for now. 
                        * TODO: Research meaning of these parameters
                        */
                        ScorerConfig config = new ScorerConfig ();
                        config.Algorithm = RecursionAlgo.Prob;
                        var diag_cross = 4;
                        var scoreDiff = 12;
                        var fastScoreThreshold = -12.5;
                        var bo = new BandingOptions (diag_cross, scoreDiff);
                        var qc = new QuiverConfig(ctx_params, bo, fastScoreThreshold);
                        config.Qconfig = qc;
                       //perf.Time(zmw.HoleNumber, "CCSMutationTesting");
                        using (var scorer = new MultiReadMutationScorer (regions, bases, initialTpl, config))//, snr))
                        {
                            result = FindConsensus.MultiReadConsensusAndQv (scorer, regions, bases.Zmw,
                                maxIterations, out iterationsTaken);
                        }

                    }
                }
            }


            #if FALSE
            // Recalibrate the accuracy prediction.  The model based estimate is not quite right.
            // We should revisit this later when we fix the model.
            float recalibratedPredAcc = RecalibratedAccuracyPrediction(bases.Zmw.Movie.SequencingChemistry, result.PredictedAccuracy, result.NumPasses);
            ((ZmwMetricsBases)result.Metrics).ReadScore = recalibratedPredAcc;
            #endif

            // Check that the read appears to be reasonably accurate based on the CCS QVs
            var passFilter = QualityFilter(result, bases);

            if ( (passFilter == CCSResultType.Success || passFilter == CCSResultType.PostCCSAccuracy) &&
                result.Sequence.Length >= Config.MinLength &&
                result.Sequence.Length <= Config.MaxLength)
            {
                // Valid CCS read
                return new Tuple<CCSResultType, IZmwConsensusBases>(CCSResultType.Success, result);
            }
            else
            {
                // We hads a bad CCS result - just bail and return the metrics of the best subread.
                return new Tuple<CCSResultType, IZmwConsensusBases>(passFilter, ZmwConsensusBases.Null(bases.Zmw));
            }
        }

        /// <summary>
        /// A filter to clean up low-quality CCS reads
        /// </summary>
        /// <returns></returns>
        public CCSResultType QualityFilter(IZmwConsensusBases consensusBases, IZmwBases bases)
        {
            // Don't emit very short CCS reads
            CCSResultType result =  consensusBases.NumBases < 5 ? CCSResultType.PostCCSShort : CCSResultType.Success ;

            var ccsAccPred = consensusBases.PredictedAccuracy;

            if (ccsAccPred < Config.MinPredictedAccuracy)
                result = CCSResultType.PostCCSAccuracy;

            // Don't emit low SNR reads
            if (bases.Metrics.HQRegionSNR.Min () < 3f)
                result = CCSResultType.OutsideSNR;
            
            // Don't emit palindromic CCS sequences
            var ccsSeq = consensusBases.Sequence;
            
            // only check for palindrome on multipass reads below a certain length
            if (consensusBases.NumPasses >= 2 && ccsSeq.Length < 4000)
            {
                // Use lo-mem global alignment -- this can cause OO
                var palindromeAl = LinearMemAlign.Align(ccsSeq, DNA.ReverseComplement(ccsSeq), GlobalSettings.FullGlobal);

                if (palindromeAl.Accuracy > 0.85) {
                    result = CCSResultType.PostCCSPalindrome;
                    Console.WriteLine ("Palindrome at hole: " + bases.Zmw.HoleNumber);
                }
            }

            // You made it!
            return result;
        }

        /// <summary>
        /// Compute SMC for one trace.
        /// </summary>
        public Tuple<TrialTemplate, float, List<AlignedSequenceReg>> GetPoaAndRegions(IZmwBases bases)
        {
            // Don't even bother if this is not a productive ZMW. 
            // This will filter out low SNR or otherwise crappy, which take a lot of time
            if(bases.Metrics.Productivity != ProductivityClass.Productive)
                return new Tuple<TrialTemplate, float, List<AlignedSequenceReg>>(null, 0.0f, new List<AlignedSequenceReg>());

            //perf.Time(zmw.HoleNumber, "CCSPartition");

            // Find adapters
            var readRegions = partitioner.GetPartition(bases, AdapterPadBases);

            // Report some stats
            //perf.Measure(zmw.HoleNumber, "TotalBases", (int)bases.NumBases);
            //perf.Measure(zmw.HoleNumber, "RawPasses", readRegions.InsertRegions.Count);
            //perf.Measure(zmw.HoleNumber, "ReadScore", bases.Metrics.ReadScore);

            // 0 inserts -- empty result
            if (readRegions.InsertRegions.Count == 0)
                return new Tuple<TrialTemplate, float, List<AlignedSequenceReg>>(null, 0.0f, new List<AlignedSequenceReg>());

            // In 2.0 we only emit >0-length sequences for true 'CCS' reads -- best estimate reads & CCS are moving to secondary
            // However, we supply the metrics for 'Read of Insert' reporting here, via the metrics attached to IZmwConsensusBases.
            // If there is not 1 full pass, just report the stats, not the sequence of the longest subread.
            // We will run CCS on everything else.  If it dones't make the 90% threshold, then we will not report any sequence for it, but we will report the metics.

            if (readRegions.InsertRegions.Count < Config.MinFullPasses)
            {
                ZmwConsensusBases.Null(bases.Zmw);
            }

            // Compute insert length
            var insertLength = readRegions.InsertRegions.Average(r => r.Length);
            //perf.Measure(zmw.HoleNumber, "AverageInsert", insertLength);

            // If this is an adapter dimer, or a very short insert then bail
            if (insertLength < 10)
                return new Tuple<TrialTemplate, float, List<AlignedSequenceReg>>(null, 0.0f, new List<AlignedSequenceReg>());


            //perf.Time(zmw.HoleNumber, "ComputePOA");

            // Quality metrics of the POA consensus
            float graphScore = 0.0f;
            //List<MutationScore> graphMutations = null;

            // The POA step must emit an initial template, and the subread regions to use
            TrialTemplate initialTpl = null;
            List<AlignedSequenceReg> regions = null;
            //List<MutationScore> startMutations = null;

            //var numCompletePasses = readRegions.InsertRegions.Count(r => r.AdapterHitBefore && r.AdapterHitAfter);

            // New-style POA setup -- use in all cases -- do a local alignment between subreads to find
            // overlapping sections
            initialTpl = FindConsensus.InitialConsensusTemplate(readRegions, bases, out graphScore,
                out regions, AdapterPadBases, Adapter);

            return new Tuple<TrialTemplate, float, List<AlignedSequenceReg>>(initialTpl, graphScore, regions);
        }


        /// <summary>
        /// Indicates that the pipeline has encountered an unrecoverable error and will shut down.
        /// Any clean-up tasks should be executed here.
        /// </summary>
        /// <param name="e"></param>
        public override void OnError(Exception e)
        {
            // No cleanup to do.
        }
    }
}
