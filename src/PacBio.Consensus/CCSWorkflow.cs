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

namespace PacBio.Consensus
{


     public class ConsensusConfig
    {
        /// <summary>
        /// Quiver parameters
        /// </summary>
        public ScorerConfig ScConfig;

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
    public class CCSStream : PipelineMapper<ConsensusConfig, IZmwBases, IZmwConsensusBases>
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

        public ScorerConfig ScConfig;


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
            ScConfig = config.ScConfig;
        }



        public override bool Map(IZmwBases bases)
        {
            var zmw = bases.Zmw;
           try
            {
                return InnerMap(bases);
            }
            catch (Exception e)
            {
                return false
            }
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

             /// <summary>
        /// Compute SMC for one trace.
        /// </summary>
        public bool InnerMap(IZmwBases bases)
        {
            // Don't even bother if this is not a productive ZMW. 
            // This will filter out low SNR or otherwise crappy, which take a lot of time
            if (bases.Metrics.Productivity != ProductivityClass.Productive)
                return false;


            if (!Config.SnrCut.Contains(bases.Metrics.HQRegionSNR))
                return false;


            // Find adapters
            var readRegions = partitioner.GetPartition(bases, AdapterPadBases);

            // 0 inserts -- empty result
            if (readRegions.InsertRegions.Count == 0)
                return false;

            // In 2.0 we only emit >0-length sequences for true 'CCS' reads -- best estimate reads & CCS are moving to secondary
            // However, we supply the metrics for 'Read of Insert' reporting here, via the metrics attached to IZmwConsensusBases.
            // If there is not 1 full pass, just report the stats, not the sequence of the longest subread.
            // We will run CCS on everything else.  If it dones't make the 90% threshold, then we will not report any sequence for it, but we will report the metics.

            if (readRegions.InsertRegions.Count(r => r.AdapterHitAfter && r.AdapterHitBefore) < Config.MinFullPasses)
            {
                return false;
            }

            // Compute insert length
            var insertLength = readRegions.InsertRegions.Average(r => r.Length);
            //perf.Measure(zmw.HoleNumber, "AverageInsert", insertLength);

            // If this is an adapter dimer, or a very short insert then bail
            if (insertLength < 10)
                return false;


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

            if (initialTpl.Length > 5) {
                return true;
            }
            return false;


        }

        /// <summary>
    
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


    }
}
