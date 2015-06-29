using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.IO.Compression;
using System.Linq;
using System.Runtime.Serialization.Formatters.Binary;
using ConsensusCore;
using PacBio.Align;
using PacBio.IO;
using PacBio.Utils;

namespace PacBio.Consensus
{

   
    /// <summary>
    /// Utility code for getting data into ConsensusCore
    /// </summary>
    public class ConsensusCoreWrap
    {


     
        private static T[] FilterArray<T>(IEnumerable<T> array, IList<bool> filter)
        {
            return array.Where((v, idx) => filter[idx]).ToArray();
        }


        private static FloatArray MakeFloatArray(IList<float> a, int start = 0, int length = -1)
        {
            if (length < 0)
                length = a.Count;

            var fa = new FloatArray(length);

            for (var i = 0; i < length; i++)
                fa.setitem(i, a[start + i]);

            return fa;

        }


        private static FloatArray MakeFloatArray(IList<byte> a, int start = 0, int length = -1)
        {
            if (length < 0)
                length = a.Count;

            var fa = new FloatArray(length);
            for (int i = 0; i < length; i++)
                fa.setitem(i, a[start + i]);

            return fa;
        }


        private static FloatArray MakeFloatArray(IList<char> a, int start = 0, int length = -1)
        {
            if (length < 0)
                length = a.Count;

            var fa = new FloatArray(length);

            for (int i = 0; i < length; i++)
                fa.setitem(i, a[start + i]);

            return fa;
        }


        private static FloatArray MakeFloatArray(IList<sbyte> a, int start = 0, int length = -1)
        {
            if (length < 0)
                length = a.Count;

            var fa = new FloatArray(length);
            for (var i = 0; i < length; i++)
                fa.setitem(i, a[start + i]);

            return fa;
        }


      
       
    }

    public class ScorerConfig : IDisposable
    {
        public QuiverConfig Qconfig;

        public RecursionAlgo Algorithm = RecursionAlgo.Prob;


        public void Dispose()
        {
            if (Qconfig != null) {
                Qconfig.Dispose ();
                Qconfig = null;
            }

        }
    }


    /// <summary>
    /// Computes the scores of a mutations with respect to an initial template for a set of reads.  This is essentially a wrapper around the MultiReadMutationScorer
    /// class in ConsensusCore.  Some things, notably the Mutation type and the 'TrialTemplate' type have one type on the C# side and a different type in ConsensusCore/SWIG.
    /// There may be some opportunity to clean that up.
    /// </summary>
    public class MultiReadMutationScorer : IDisposable
    {
        private static PacBioLogger logger = PacBioLogger.GetLogger("MultiReadMutationScorer");

        private static void Log(LogLevel level, string msg)
        {
            logger.Log(level, msg);
        }

        private static void Log(LogLevel level, string msg, params object[] args)
        {
            logger.Log(level, String.Format(msg, args));
        }

        internal static ByteVector MakePWVector(IList<ushort> a, int start, int length = -1)
        {
            if (length < 0)
                length = a.Count;
            /* Despite what the SWIG documentation might indicate, under the hood
               new ByteVector(length) performs the following (setting capacity and not size):
                        pv = new std::vector< unsigned char >();
                        pv->reserve(length);
                and so attempts to index it directly via [] leads to an out of bounds exception.
                Will use .Add instead for now.
            */
            var ba = new ByteVector (length);
            int i = 0;              
            for (i = 0; i < length; i++) {
                var value = a [start + i];
                if (value <= 20)
                    ba.Add((byte)value);
                else
                    ba.Add((byte)20);
            }           
            return ba;
        }
        internal static ByteVector MakeMqvVector(IZmwBases bases, out string seq_a, int start, int length = -1)
        {
            List<char> seq = new List<char> (length);
            List<byte> mqv = new List<byte> (length);

            var delTags = bases.DeletionTag.ToArray ();
            var mqvs = bases.MergeQV.ToArray ();
            var oseq = bases.Sequence.Substring (start, length);
                         
            for (int i = 0; i < length; i++) {
                if (i < (length - 1)) {
                    var delTag = delTags [start + i];
                    if (delTag != 'N' && delTag != '-') {
                        seq.Add (delTag);
                        mqv.Add ((byte)0);
                    }
                }
                var value = mqvs [start + i];
                seq.Add (oseq [i]);
                var val2 = Convert.ToInt32(Math.Floor(((double) value) / 3.0)) + 1;
                if (val2 <= 20)
                    mqv.Add((byte)val2);
                else
                    mqv.Add((byte)20);
            }

            var mqv_a = new ByteVector (mqv.Count);
            foreach(byte b in mqv) {
                mqv_a.Add(b);
            }
            seq_a = new String (seq.ToArray ());
            return mqv_a;
        }

        internal static ByteVector MakeIqvVector(IList<byte> a, int start, int length = -1)
        {
            if (length < 0)
                length = a.Count;
            /* Despite what the SWIG documentation might indicate, under the hood
               new ByteVector(length) performs the following (setting capacity and not size):

                        pv = new std::vector< unsigned char >();
                        pv->reserve(length);
                        
                and so attempts to index it directly via [] leads to an out of bounds exception.
                Will use .Add instead for now.
            */
            var ba = new ByteVector (length);
            int i = 0;              
            for (i = 0; i < length; i++) {
                var value = a [start + i];
                var val2 = Convert.ToInt32(Math.Floor(((double) value) / 3.0)) + 1;
                if (val2 <= 20)
                    ba.Add((byte)val2);
                else
                    ba.Add((byte)20);
            }           
            return ba;
        }

        public float ScoreDiff = 15;
        public float DynamicAdjustFactor = 0.0f;
        public float DynamicAdjustOffset = -0.3f;
        #if DIAGNOSTIC
        public const float AddThreshold = 0.15f;
        #else
        public const float AddThreshold = 1.0f;
        #endif

        /// <summary>
        /// Number of bases of the adapter included at front of sequence
        /// </summary>
        private int StartAdapterBases;

        /// <summary>
        /// Number of bases of adapter included at the end of sequence
        /// </summary>
        private int EndAdapterBases;

        public SparseQvSumProductMultiReadMutationScorer scorer;


        public List<int[]> delTags;
        /// <summary>
        /// Construct a mutation evalutor for the pulse observations in encapsulated by read, on template TrialTemplate.
        /// </summary>
        public MultiReadMutationScorer(IEnumerable<AlignedSequenceReg> regions, IZmwBases bases,
            TrialTemplate trialTemplate, ScorerConfig config)
        {
            StartAdapterBases = trialTemplate.StartAdapterBases;
            EndAdapterBases = trialTemplate.EndAdapterBases;

            var strandTpl = trialTemplate.GetSequence(Strand.Forward);

            if (config.Algorithm == RecursionAlgo.Viterbi)
            {
                throw new Exception ("Integrate don't optimize");
                //scorer = new SparseQvViterbiMultiReadMutationScorer(config.Qconfig, strandTpl);
            }
            else if (config.Algorithm == RecursionAlgo.Prob)
            {
                scorer = new SparseQvSumProductMultiReadMutationScorer(config.Qconfig, strandTpl);
            }
            else
            {
                throw new ApplicationException("unrecognized recursion algorthim");
            }

//            StreamWriter sw = new StreamWriter ("Output.cpp");
//            sw.WriteLine ("SNR snr(" + snr.A.ToString () + "," + snr.C + "," + snr.G + "," + snr.T + ");");
//            sw.WriteLine ("ContextParameters ctx_params(snr);");
//            sw.WriteLine ("QuiverConfig qc(ctx_params, bo, fast_sçore_threshold);\nSparseSimpleSumProductMultiReadMutationScorer scorer(qc, temp);");
//
            delTags = new List<int[]>();
            foreach (var r in regions)
            {
                var name = TraceReference.CreateSpringfieldSubread(bases, r.Start, r.End);
                //var seq = bases.Sequence.Substring (r.Start, r.End - r.Start);
                //using (var qsf = ConsensusCoreWrap.MakeQvSequenceFeatures(r.Start, r.End - r.Start, bases))

                string seq;
                using (var iqvs = MakeMqvVector(bases, out seq, r.Start, r.End - r.Start))  
                using (var read = new Read(name, seq, iqvs, iqvs)) // Hack to make sure we use merge QV
                using (var mappedRead = new MappedRead(read, (StrandEnum) r.Strand, r.TemplateStart, r.TemplateEnd,
                                                       r.AdapterHitBefore, r.AdapterHitAfter))
                {
                    var result = scorer.AddRead (mappedRead, AddThreshold);
                    //Console.WriteLine (scorer.BaselineScore ());
                    if (result != AddReadResult.SUCCESS) {
                        
                        
                        #if DIAGNOSTIC
                        scorer.AddRead(mappedRead, 1.0f);

                        var prefix = name.Replace("/", ".");

                        WriteAlphaBetaMatrices(readBitmap.Count, prefix);
                        #else
                        Log (LogLevel.DEBUG, "Skipped adding read '{0}' -- more than {1}% of entries used.", name, AddThreshold * 100);
                        #endif
                    } else {
                        delTags.Add(bases.DeletionTag.Skip (r.Start).Take (r.End - r.Start).Select(z=>(int)z).ToArray());

                    }
                }
            }
        }

        int ReadCounter =0;
        void DumpData(SparseQvSumProductMultiReadMutationScorer scorer, StreamWriter sw)
        {
            ReadCounter++;
            var cnt = scorer.BaselineScores ().Count;
            for (int i = 0; i < cnt; i++) {
                var rnd = scorer.GetRead (i);
                if (rnd.IsActive) {
                    var rname = "r" + i;
                    string tmpName = "temp_" + rname;
                    sw.WriteLine ("std::string " + tmpName + " = \"" + rnd.Scorer.Template().GetTemplateSequence() + "\";");

                    var rl = "Read " + rname + "(\"" + rnd.Read.Name + "\",\"" + rnd.Read.Sequence + "\");";
                    sw.WriteLine (rl);
                    var mrname = "mr" + ReadCounter;
                    var strand = "StrandEnum::" + (rnd.Read.Strand == StrandEnum.REVERSE_STRAND ? "REVERSE_STRAND" : "FORWARD_STRAND");
                    var mrl = "MappedRead " + mrname + "(r" + ReadCounter + ", " +
                        strand + ", " + rnd.Read.TemplateStart+ ", " + rnd.Read.TemplateEnd + " );";
                    sw.WriteLine (mrl);

                    var next = "auto result" + ReadCounter + " = scorer.AddRead(" + mrname + ", 1.0);";
                    sw.WriteLine (next);

                    var iqv = "std::vector<unsigned char> iqvs = {" + String.Join (",", rnd.Read.Iqvs.Select (x => x.ToString ())) + "};";
                    sw.WriteLine (iqv);

                    sw.Flush ();
                }
            }


        }

        void DumpData(Read rd, MappedRead mr, StrandEnum se, int ts, int tend, bool hitbefore, bool hitafter, StreamWriter sw, string tpl)
        {
            ReadCounter++;
            sw.WriteLine("std::string temp = \"" + tpl.Substring(mr.TemplateStart, mr.TemplateEnd - mr.TemplateStart +1) + "\";" );

            var rname = "r" + ReadCounter;
            var rl = "Read "+rname + "(\"" + rd.Name + "\",\"" + rd.Sequence + "\");";
            sw.WriteLine (rl);
            var mrname = "mr" + ReadCounter;
            var strand = "StrandEnum::" + (se == StrandEnum.REVERSE_STRAND ? "REVERSE_STRAND" : "FORWARD_STRAND");
            var mrl = "MappedRead " + mrname + "(r" + ReadCounter + ", " +
                strand + ", " + ts + ", " + tend + ", " + hitbefore.ToString().ToLower() + ", " + hitafter.ToString().ToLower() + " );";
            sw.WriteLine (mrl);

            var next = "auto result" + ReadCounter + " = scorer.AddRead(" + mrname + ", 1.0);";
            sw.WriteLine (next);

            var iqv = "std::vector<unsigned char> iqvs = {" + String.Join (",", mr.Iqvs.Select (x => x.ToString ())) + "};";
            sw.WriteLine (iqv);

            sw.Flush ();


        }




        /// <summary>
        /// Gets the number of allocated cells for each added MappedRead
        /// </summary>
        public int[] AllocatedEntries
        {
            get
            {
                using (var allocatedMatrixEntries = scorer.AllocatedMatrixEntries())
                {
                    var r = new int[allocatedMatrixEntries.Count];
                    allocatedMatrixEntries.CopyTo(r);
                    return r;
                }
            }
        }

        /// <summary>
        /// Gets the used cells for each added MappedRead
        /// </summary>
        public int[] UsedEntries
        {
            get
            {
                using (var usedMatrixEntries = scorer.UsedMatrixEntries())
                {
                    var r = new int[usedMatrixEntries.Count];
                    usedMatrixEntries.CopyTo(r);
                    return r;
                }
            }
        }


        private void WriteCsvMatrix(int i, bool alpha, string prefix)
        {
            using (var stream = new FileStream(prefix + ".csv.gz", FileMode.Create))
            using (var gzstream = new GZipStream(stream, CompressionMode.Compress))
            using (var writer = new StreamWriter(gzstream))
            {
                var mat = alpha ? scorer.AlphaMatrix(i) : scorer.BetaMatrix(i);
                var nrow = mat.Rows();
                var ncol = mat.Columns();

                for (int x = 0; x < nrow; x++)
                {
                    var row = Enumerable.Range(0, ncol)
                        .Select(y => mat.IsAllocated(x, y) ? mat.Get(x, y) : float.NaN);
                    writer.WriteLine(row.Select(v => v.ToString()));
                }
            }
        }


        public void WriteAlphaBetaMatrices(int i, string prefix)
        {
            WriteCsvMatrix(i, true, prefix + ".alpha");
            WriteCsvMatrix(i, false, prefix + ".beta");
        }


        /// <summary>
        /// Estimate of the memory consumption in MB of the recursion matrices
        /// </summary>
        public float AllocatedSizeMB
        {
            get
            {
                return  AllocatedEntries.Sum() * 4.0f / 1e6f;
            }
        }

        /// <summary>
        /// Number of forward-backward flip-flops taken to converge for each read
        /// </summary>
        public int[] NumFlipFlops
        {
            get
            {
                using (var nFlipFlops = scorer.NumFlipFlops())
                {
                    var r = new int[nFlipFlops.Count];
                    nFlipFlops.CopyTo(r);
                    return r;
                }
            }
        }

        /// <summary>
        /// Access the ConsensusCore MappedRead object for each read
        /// </summary>
        public MappedRead[] MappedReads
        {
            get
            {
                var r = new MappedRead[scorer.NumReads ()];
                for (int i = 0; i < scorer.NumReads(); i++)
                {
                    r[i] = scorer.Read(i);
                }

                return r;
            }
        }

        public List<ReadStateType> MappedStates {
            get {
                List<ReadStateType> toR = new List<ReadStateType> (scorer.NumReads ());
                for (int i = 0; i < scorer.NumReads(); i++) {
                    var r = scorer.GetRead (i);
                    if (r != null) {
                        toR.Add (scorer.GetRead (i));
                    }
                }
                return toR;
            }
        }

        public int NumReads
        {
            get { return scorer.NumReads(); }
        }
        #if FALSE
        /// <summary>
        /// Copy constructor
        /// </summary>
        private MultiReadMutationScorer(MultiReadMutationScorer other)
        {
            StartAdapterBases = other.StartAdapterBases;
            EndAdapterBases = other.EndAdapterBases;

            if (other.scorer is SparseSseQvMultiReadMutationScorer)
            {
                scorer = new SparseSseQvMultiReadMutationScorer((SparseSseQvMultiReadMutationScorer) other.scorer);
            }
            else if (other.scorer is SparseSseQvSumProductMultiReadMutationScorer)
            {
                scorer = new SparseSseQvSumProductMultiReadMutationScorer((SparseSseQvSumProductMultiReadMutationScorer) other.scorer);
            }
            else
            {
                throw new ApplicationException("Unrecognized MultiReadMutationScorer type");
            }
        }

        /// <summary>
        /// Make a clone of the scorer with the same set of reads.
        /// </summary>
        public MultiReadMutationScorer Clone()
        {
            // Invoke copy constructor
            return new MultiReadMutationScorer(this);
        }
        #endif

        public double BaselineScore()
        {
            return scorer.BaselineScore();
        }


        /// <summary>
        /// Return the ratio of the P(T'|O) / P(T|O) where T is the template that the MutationEvaluator was
        /// constructed with, and T' is the template resulting from the application of Mutation m to template T.
        /// </summary>
        /// <param name="m">The mutation to apply to the original template</param>
        /// <returns>The ratio of the mutated and original template likelihoods</returns>
        public double ScoreMutation(Mutation m)
        {
            return scorer.Score(m.Type, m.TemplatePosition, m.Base);
        }

        public double ScoreMutationFast(Mutation m)
        {
            using (var ccMutation = m.ToConsensusCoreType())
            {
                return scorer.FastScore(ccMutation);
            }
        }

        public double ScoreMutationWeighted(Mutation m, double[] mappingRatios)
        {
            var scores = GetScores(m);

            double acc = 0.0;

            for (int i = 0; i < scores.Length; i++)
            {
                if (mappingRatios[i] > 1e7)
                {
                    // Read maps essentially uniquely to this haplotype.
                    acc += scores[i];
                }
                else if (mappingRatios[i] < 1e-7)
                {
                    // Doesn't map strongy enough to affect outcome.
                }
                else
                {
                    // Partial mapping
                    var expScore = (Math.Exp(scores[i]) + (1.0f/mappingRatios[i])) / (1.0f + (1.0f/mappingRatios[i]));
                    var score = Math.Log(expScore);
                    acc += (double) score;
                }
            }

            return acc;
        }

        /// <summary>
        /// Score a compound mutation by packing the single base mutations into a single multi-base mutation.
        /// FIXME - not debugged in ConsensusCore yet.
        /// </summary>
        public double ScoreCompoundMutation(List<Mutation> mutations)
        {
            var ccMuts = mutations.Map(m => new ConsensusCore.Mutation(m.Type, m.TemplatePosition, m.Base));

            try
            {
                var mutVect = new MutationVector();
                foreach (var m in ccMuts)
                {
                    mutVect.Add(m);
                }

                return 0.0f;

                // FIXME support 'compound' mutations
                // var combinedMut = ConsensusCore.Mutation.Combine(mutVect, scorer.Template());
                // return scorer.Score(combinedMut);
            }
            finally
            {
                // Clean up C++ mutation objects
                ccMuts.ForEach(m => m.Dispose());
            }
        }

        /// <summary>
        /// Get the Quiver score deltas for this mutation for each read with a default value for unobserved mutations
        /// </summary>
        public double[] GetScoresWithDefault(Mutation m, double defaultValue)
        {
            using (var scores = scorer.Scores(m.Type, m.TemplatePosition, m.Base, defaultValue))
            {
                var r = new double[scores.Count];
                scores.CopyTo(r);
                return r;
            }
        }

        /// <summary>
        /// Aggregate the quiver score deltas according to some supplied function
        /// </summary>
        public double AggregateScoresWithDefault(Mutation m, double defaultValue, Func<double, double, double> f)
        {
            using (var scores = scorer.Scores(m.Type, m.TemplatePosition, m.Base, defaultValue))
            {
                return scores.Aggregate(defaultValue, f);
            }
        }

        /// <summary>
        /// Get the Quiver score deltas for this mutation for each read.
        /// </summary>
        public double[] GetScores(Mutation m)
        {
            return GetScoresWithDefault(m, 0);
        }

        /// <summary>
        /// Get the baseline Quiver score of each read.
        /// </summary>
        public double[] GetBaselineScores()
        {
            using (var scores = scorer.BaselineScores())
            {
                var r = new double[scores.Count];
                scores.CopyTo(r);
                return r;
            }
        }

        /// <summary>
        /// Apply a set of mutations to the template
        /// </summary>
        public void ApplyMutations(List<Mutation> mutations)
        {
            // Convert to the CC mutation type
            // NOTE - you must ensure that the CC mutations are not GC'd/deleted before the C++ ApplyMutations method is finished.
            var ccMuts =
                mutations.Map(m => new ConsensusCore.Mutation(m.Type, m.TemplatePosition, m.Base));
            try
            {
                ApplyToScorer(ccMuts);
            }
            catch (Exception e)
            {
                Console.WriteLine(e.Message);
                if(e.InnerException != null)
                    Console.WriteLine(e.InnerException.Message);
                Console.WriteLine("Got bad mutations!");
                //System.Diagnostics.Debugger.Break();
            }

            // Clean up C++ mutation objects
            ccMuts.ForEach(m => m.Dispose());
        }

        private void ApplyToScorer(IEnumerable<ConsensusCore.Mutation> muts)
        {
            var mutVect = new MutationVector();
            foreach(var m in muts)
                mutVect.Add(m);

            scorer.ApplyMutations(mutVect);
        }

        /// <summary>
        /// Current template
        /// </summary>
        public TrialTemplate Template
        {
            get { 

                return new TrialTemplate(scorer.Template (StrandEnum.FORWARD_STRAND).tpl, StartAdapterBases, EndAdapterBases); 
                //throw new Exception();// 
            }
               
        }


        /// <summary>
        /// Performs application-defined tasks associated with freeing, releasing, or resetting unmanaged resources.
        /// </summary>
        public void Dispose()
        {
            if (scorer != null)
            {
                scorer.Dispose();
                scorer = null;
            }
        }

        /// <summary>
        /// Static helper method for writing out a Quiver score matrix.  Useful for making test datasets for phasing.
        /// </summary>
        public static void SaveMultiSplitScores(MultiReadMutationScorer sc, string filename)
        {
            // Use all possible substitutions as potential splitting mutations
            var muts =
                GenerateMutations.GenerateUniqueMutations(sc.Template)
                    .Where(s => s.Type == MutationType.SUBSTITUTION)
                    .ToArray();

            var mutScoreVectors = muts.Map(sc.GetScores);
            var positions = muts.Map(m => m.TemplatePosition);

            var binaryFormatter = new BinaryFormatter();
            using (
                var serializationStream = new FileStream(filename, FileMode.Create,
                    FileAccess.Write,
                    FileShare.None))
            {
                binaryFormatter.Serialize(serializationStream, new Tuple<double[][], int[]>(mutScoreVectors, positions));
            }
        }
    }
}


