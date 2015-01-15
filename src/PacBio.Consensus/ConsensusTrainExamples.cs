using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;
using PacBio.IO;
using PacBio.Utils;
using PacBio.Align;

namespace PacBio.Consensus
{
    /// <summary>
    /// Associate cmp.h5 alignment data with bas.h5 read data
    /// </summary>
	public class TraceSet
	{
		OldCmpH5Reader cmp;
		BasCollection bas;
        SnrCut snrCut;

        public TraceSet(string cmpFile, string fofn, SnrCut snrCut = null)
		{
            cmp = OldCmpH5Reader.CreateReader(cmpFile);
			bas = BasCollection.FromFofn(fofn);
            this.snrCut = snrCut ?? SnrCut.PassAll;
		}

		public IEnumerable<Trace> Traces 
		{
			get
			{
                return bas.Reads
                    .Where(r => snrCut.Contains(r.Metrics.HQRegionSNR))
                    .Select (MakeTrace);
			}
		}

		public int NumTraces
		{
			get {
                return Traces.Count();
			}
		}

		public Trace MakeTrace(IZmwBases b)
		{
            var alns = cmp.GetAlignmentsForZmw(b.Zmw.Movie.MovieName, b.Zmw.HoleNumber).ToArray();

            IAlnSummary bestAln = null;

            if (alns.Length > 0)
                bestAln = alns.OrderByDescending(a => Math.Pow (a.Accuracy, 2.0) * a.TemplateLength).First();

			return new Trace {
				ZmwBases = b,
				MultiAlignment = alns,
				SmithWatermanAlignment = bestAln,
                ReportsFolder = cmp.ReportsFolder
			};
		}
	}

	public class Trace
	{
        public OldCmpH5Reader Reader { get; set; }
        public string ReportsFolder { get; set; }

		public IZmwBases ZmwBases { get; set; }
        public IAlnSummary[] MultiAlignment { get; set; }
        public IAlnSummary SmithWatermanAlignment { get; set; }

	    public int HoleNumber
	    {
            get { return ZmwBases.Zmw.HoleNumber; }
	    }

	    public IMovieMetadata Metadata
	    {
            get { return ZmwBases.Zmw.Movie;  }
	    }

        public IAlignment GetAlignment(IAlnSummary s)
        {
            return Reader.ReadAlignment(s);
        }
	}

    public class CCSExample
    {
        public Trace Trace;
        public string Reference;
        public TrialTemplate CorrectTrialTemplate;
        public AlignedSequenceReg[] Regions;

        public static CCSStream ccsAlgo = CCSStream.DefaultConfig;

        public static SimpleAlignment SemiGlobalAlign(string refSeg, string read)
        {
            var penalty = new short[] { -7, -7, -13, 4, -4 };
            //var refSeg = reference.Substring(refStart, refEnd - refStart + 1);
            var a1 = GlobalAlign.Align(refSeg, read, penalty, false, true);

            var rc = DNA.ReverseComplement(refSeg);
            var a2 = GlobalAlign.Align(rc, read, penalty, false, true, Strand.Reverse);
            if (a1.Accuracy > a2.Accuracy)
            {
                return a1;
            }
            
            return a2;
        }

        public static SimpleAlignment LocalAlign(string refSeg, string read)
        {
            var penalty = new short[] { -7, -7, -13, 3, -4 };
            var a1 = CircularSmithWaterman.Align(refSeg, read, penalty);

            var rc = DNA.ReverseComplement(refSeg);
            var a2 = CircularSmithWaterman.Align(rc, read, penalty);
            if (a1.TemplateLength > a2.TemplateLength)
            {
                return a1;
            }

            return a2;

        }

        public static SimpleAlignment FullGlobalAlign(string refSeg, string read)
        {
            var penalty = new short[] { -7, -7, -13, 4, -4 };
            //var refSeg = reference.Substring(refStart, refEnd - refStart + 1);
            var a1 = GlobalAlign.Align(refSeg, read, penalty, true, true);

            var rc = DNA.ReverseComplement(refSeg);
            var a2 = GlobalAlign.Align(rc, read, penalty, true, true, Strand.Reverse);

            if (a1.Accuracy > a2.Accuracy)
            {
                return a1;
            }

            return a2;
        }


        public enum ExampleMode
        {
            Sample,
            Fast
        }

        /// <summary>
        /// Gets a simple set of examples, used for comparing parameters.
        /// </summary>
        /// <returns>The simple examples.</returns>
        /// <param name="traceSet">Trace set.</param>
        /// <param name="referenceContigs">Reference contigs.</param>
        /// <param name="examplesToGet">Examples to get.</param>
        /// <param name="curReference">Current reference.</param>
        public static List<CCSExample> GetSimpleExamples(TraceSet traceSet, Dictionary<string, string> referenceContigs, 
            int examplesToGet, string curReference)
        {           
            int totalNeeded = examplesToGet;
            var toReturn = new List<CCSExample> (totalNeeded);
           
            var ccsTraces = traceSet.Traces.Where (z => z.SmithWatermanAlignment!= null && z.SmithWatermanAlignment.ReferenceName == curReference );
            var nTotal = traceSet.NumTraces;
            var nTried = 0;
            int accepted = 0;
            var rejects = new Dictionary<string, int>
            {
                {"UnAligned", 0}, {"Al70Acc80", 0}, {"ETControl", 0},
                {"Passes<03", 0}, {"Passes>80", 0}, {"BadSnrChk", 0},
                {"WeirdAlgn", 0}, {"NotOkay",0}
            };

            foreach(var t in ccsTraces) 
            {
                nTried++;
                if(t.MultiAlignment.Length < 2) {
                    ++rejects["UnAligned"];
                    continue;
                }
                // Cap the accuracy at 80% - so we will get the longest read among those above 80% accuracy.
                var bestAl = t.MultiAlignment.OrderByDescending
                    (al => {
                        var accCap = Math.Min(al.Accuracy, 0.85);
                        return accCap * al.TemplateLength;
                    }).First();

                if(bestAl.TemplateLength < 70 || bestAl.Accuracy < 0.80) {
                    ++rejects["Al70Acc80"];
                    continue;
                }
                if(bestAl.ReferenceName.Contains("ET")) {
                    ++rejects["ETControl"];
                    continue;
                }

                if(!CheckSnr(t)) {
                    ++rejects["BadSnrChk"];
                    continue;
                }
                var r = ccsAlgo.GetPoaAndRegions(t.ZmwBases);
                var passes = r.Item3;
                var tpl = r.Item1;
                var poaScore = r.Item2;

                if(r.Item1 == null)
                {
                    continue;
                }
                // Want to work on at least 3 passes
                if(passes.Count < 3) {
                    ++rejects["Passes<03"];
                    continue;
                }
                // Don't use more than 80 passes -- a waste of time because the error rate should be low
                if(passes.Count > 80) {
                    ++rejects["Passes>80"];
                    continue;
                }

                var refStart = bestAl.TemplateStart;
                var refLength = bestAl.TemplateLength;

                var refSeq = referenceContigs[bestAl.ReferenceName];

                const int templateExtend = 100;
                var start = Math.Max(refStart - templateExtend, 0);
                var len = Math.Min(refSeq.Length - 1 - start, refLength + 2 * templateExtend);

                var rref = refSeq.Substring(start, len);
                var poaAl = CCSExample.SemiGlobalAlign(rref, tpl.Sequence);

                var badAlign = MultiAlignCheck(t);

                if(poaAl.Accuracy < 0.87) {
                    var err = badAlign.Item1;
                    var msg = badAlign.Item2;
                    Console.WriteLine(@"Skipping trace. POA Acc: {0}. POAScore: {1}. Had weird alignments: {2}",
                        poaAl.Accuracy, poaScore, err ? msg : "False");
                    ++rejects["WeirdAlgn"];
                    continue;
                }

                // Construct the correct trial template!
                var rcRef = DNA.ReverseComplement(rref);
                var refChunk = poaAl.Strand == Strand.Forward ? rref : rcRef;
                var correctInsertSequence = refChunk.Substring(poaAl.TemplateStartBase, poaAl.TemplateLength);
                var rcCorrectInsert = DNA.ReverseComplement(correctInsertSequence);

                var trialTemplate = new TrialTemplate() {
                    Sequence = correctInsertSequence,
                    StartAdapterBases = 5,
                    EndAdapterBases = 5
                };

                // Remap regions to this sequence
                var newRegions = passes.Map(p => MapRegionToRef(p, t.ZmwBases, correctInsertSequence, rcCorrectInsert));

                Console.WriteLine(@"Accepting trace. POA Acc: {0}. POA Score: {1}", poaAl.Accuracy, poaScore);
                ++accepted;
                var example = new CCSExample {
                    Trace = t,
                    Reference = rref,
                    CorrectTrialTemplate = trialTemplate,
                    Regions = newRegions
                };
                toReturn.Add (example);
                if (accepted >= totalNeeded) {
                    break;
                } 
            }

            int nRejects = rejects.Sum(v => v.Value);

            Console.WriteLine(@"Reject Counts:");
            rejects.ForEach(r => Console.WriteLine(@"{0} = {1}", r.Key, r.Value));

            Console.WriteLine(@"Total: {0}", nTotal);
            Console.WriteLine(@"Viewed: Accepted[{0}] + Rejected[{1}] = {2}",
                accepted, nRejects, accepted + nRejects);
           return toReturn;
        }
        public static TrainingDataStore GetExamples(TraceSet traceSet, Dictionary<string, string> referenceContigs, 
                                                    int examplesPerReference, ReadConfigurationAssigner rca, ExampleMode mode = ExampleMode.Sample)
        {           
            //Console.WriteLine("Checking references");
			//ScanSets.VerifyReferences(rawCmpH5);

            //Console.WriteLine("Loading cmp.h5:");
			//varscans = ScanSets.FromCmpH5(rawCmpH5);

            var tds = new TrainingDataStore (rca, examplesPerReference);
            int totalNeeded = 2 * referenceContigs.Count * examplesPerReference;
            HashSet<string> okayRefs = new HashSet<string>(referenceContigs.Keys);

            //var ccsTraces = scans.SelectMany(s => s.LazyTraces).Where(CCSTraceFilter);
            var ccsTraces = traceSet.Traces.Where (z => z.SmithWatermanAlignment!= null && okayRefs.Contains(z.SmithWatermanAlignment.ReferenceName));
			var nTotal = traceSet.NumTraces;
            var nTried = 0;
            int accepted = 0;
            var rejects = new Dictionary<string, int>
                              {
                                  {"UnAligned", 0}, {"Al70Acc80", 0}, {"ETControl", 0},
                                  {"Passes<03", 0}, {"Passes>80", 0}, {"BadSnrChk", 0},
                {"WeirdAlgn", 0}, {"NotOkay",0}, {"TooManyRegions",0}
                              };

            var acceptedAsExamples = 0;
            Parallel.ForEach(ccsTraces, t =>  
            {
                try {
                    if (acceptedAsExamples >= totalNeeded) {
                        return;
                    } 
                    nTried++;

                    if(t.MultiAlignment.Length < 2) {
                        ++rejects["UnAligned"];
                    return;
                    }
                    if (t.MultiAlignment.Length > 25) {
                        ++rejects ["TooManyRegions"];
                        return;
                    }

                    // Cap the accuracy at 80% - so we will get the longest read among those above 80% accuracy.
                    var bestAl = t.MultiAlignment.OrderByDescending
                            (al => {
                        var accCap = Math.Min(al.Accuracy, 0.85);
                        return accCap * al.TemplateLength;
                    }).First();

                    if(bestAl.TemplateLength < 70 || bestAl.Accuracy < 0.80) {
                        ++rejects["Al70Acc80"];
                    return;
                    }
                    if(bestAl.ReferenceName.Contains("ET")) {
                        ++rejects["ETControl"];
                        return;
                    }

                    if(!CheckSnr(t)) {
                        ++rejects["BadSnrChk"];
                        return;
                    }
                    var r = ccsAlgo.GetPoaAndRegions(t.ZmwBases);
                    var passes = r.Item3;
                    var tpl = r.Item1;
                    var poaScore = r.Item2;

                    if(r.Item1 == null)
                    {
                        return;
                    }
                    // Want to work on at least 3 passes
                    if(passes.Count < 3) {
                        ++rejects["Passes<03"];
                        return;
                    }
                    // Don't use more than 80 passes -- a waste of time because the error rate should be low
                    if(passes.Count > 80) {
                        ++rejects["Passes>80"];
                        return;
                    }

                    var refStart = bestAl.TemplateStart;
                    var refLength = bestAl.TemplateLength;

                    var refSeq = referenceContigs[bestAl.ReferenceName];

                    const int templateExtend = 100;
                    var start = Math.Max(refStart - templateExtend, 0);
                    var len = Math.Min(refSeq.Length - 1 - start, refLength + 2 * templateExtend);

                    var rref = refSeq.Substring(start, len);
                    var poaAl = CCSExample.SemiGlobalAlign(rref, tpl.Sequence);

                    var badAlign = MultiAlignCheck(t);

                    if(poaAl.Accuracy < 0.87) {
                        var err = badAlign.Item1;
                        var msg = badAlign.Item2;
                        Console.WriteLine(@"Skipping trace. POA Acc: {0}. POAScore: {1}. Had weird alignments: {2}",
                            poaAl.Accuracy, poaScore, err ? msg : "False");
                        ++rejects["WeirdAlgn"];
                        return;
                    }

                    // Construct the correct trial template!
                    var rcRef = DNA.ReverseComplement(rref);
                    var refChunk = poaAl.Strand == Strand.Forward ? rref : rcRef;
                    var correctInsertSequence = refChunk.Substring(poaAl.TemplateStartBase, poaAl.TemplateLength);
                    var rcCorrectInsert = DNA.ReverseComplement(correctInsertSequence);

                    var trialTemplate = new TrialTemplate() {
                        Sequence = correctInsertSequence,
                        StartAdapterBases = 5,
                        EndAdapterBases = 5
                    };

                    // Remap regions to this sequence
                    var newRegions = passes.Map(p => MapRegionToRef(p, t.ZmwBases, correctInsertSequence, rcCorrectInsert));

                    Console.WriteLine(@"Accepting trace. POA Acc: {0}. POA Score: {1}. Ref: {2}", poaAl.Accuracy, poaScore, t.SmithWatermanAlignment.ReferenceName);
                    ++accepted;
                    var example = new CCSExample {
                        Trace = t,
                        Reference = rref,
                        CorrectTrialTemplate = trialTemplate,
                        Regions = newRegions
                    };
                    var success = tds.AddExample(example);
                if (success) {
                    Console.WriteLine ("Accepted");
                        Interlocked.Increment(ref acceptedAsExamples);
                } else {
                    Console.WriteLine ("Not Accepted");
                }
                    }
                    catch(Exception thrown) {
                        Console.WriteLine("ERROR\n\n\n"+thrown.Message+"\n\n\nStack Trace\n\n" + thrown.StackTrace);
                    }
               
                });

            int nRejects = rejects.Sum(v => v.Value);

            Console.WriteLine(@"Reject Counts:");
            rejects.ForEach(r => Console.WriteLine(@"{0} = {1}", r.Key, r.Value));

            Console.WriteLine(@"Total: {0}", nTotal);
            Console.WriteLine(@"Viewed: Accepted[{0}] + Rejected[{1}] = {2}",
                              accepted, nRejects, accepted + nRejects);
            Console.WriteLine ("Accepted as examples: " + acceptedAsExamples);
            return tds;
        }

        public static AlignedSequenceReg MapRegionToRef(AlignedSequenceReg reg, IZmwBases bases, string refChunk, string rcRefChunk)
        {
            var readSeq = bases.Sequence.Substring(reg.Start, reg.Length);
            var tplSeq = (reg.Strand == Strand.Forward) ? refChunk : rcRefChunk;
            //var regAl = CCSExample.SemiGlobalAlign(tplSeq, readSeq);

            var penalty = new short[] { -7, -7, -13, 4, -4 };
            var regAl = GlobalAlign.Align(tplSeq, readSeq, penalty, false, true);

            if(reg.Strand == Strand.Forward) {
                return new AlignedSequenceReg(reg.Start, reg.End, regAl.TemplateStartBase, regAl.TemplateStartBase + regAl.TemplateLength, reg.Strand);
            } else {

                var tplLen = regAl.TemplateLength;
                var tplEnd = regAl.TemplateStartBase + regAl.TemplateLength;

                return new AlignedSequenceReg(reg.Start, reg.End, tplSeq.Length - tplEnd, tplSeq.Length - regAl.TemplateStartBase, reg.Strand);
            }
        }


        public static bool CheckSnr(Trace t)
        {
            // Fast path SNR check -- see if the HQRegionSNR Metric passes -- if so then we should be GTG
			var minHqRegionSnr = t.ZmwBases.Metrics.HQRegionSNR.Min();

            if (minHqRegionSnr > 3.5)
                return true;

			return false;
        }

        public static Tuple<bool,string> MultiAlignCheck(Trace trc)
        {
            var ma = trc.MultiAlignment;

            if (ma.Length < 2)
                return new Tuple<bool, string> (false, null);

            // Check that all the alignments come from the same location
            var seedAl = ma.OrderByDescending (al => al.TemplateLength).First ();
            var maxDist = seedAl.TemplateLength * 2;

            if (ma.Any (compAl => Math.Abs (seedAl.TemplateStart - compAl.TemplateStart) > maxDist))
                return new Tuple<bool, string> (true, "Multiple reference locations");


            var hit0 = ma [0];
            var firstStrand = ma [0].Strand;
            var ada = trc.ZmwBases.AdapterHits();

            // Check that the template strand alternates properly
            for (int i = 1; i < ma.Length; i++) {
                var hitI = ma [i];

                var expectedStrand =
                    ada.Where (a =>
                        a.Start >= hit0.ReadStart + hit0.ReadLength - 1 &&
                        a.End <= hitI.ReadEnd).Count () % 2 == 0
                        ? firstStrand
                        : firstStrand.Opposite ();

                if (ma [i].Strand != expectedStrand)
                    return new Tuple<bool, string> (true, "Incorrect strand alternation");
            }

            return new Tuple<bool, string> (false, null);
        }



		public static IEnumerable<TResult> SparseTargetSample<T, TResult>(IEnumerable<T> data, int target, int total, Func<T, TResult> chooser) where  TResult : class
		{
		    return data.AsParallel().Select(chooser).Where(v => v != null).Take(target);
		    /*
			var r = new Random();

			// Number we ran the function on
			int nTested = 0;
			// Number that passed the test
			int nPassed = 0;
			// Number returned out of function
			int nReturned = 0;
			// Number of data points read from input
			int nIterated = 0;

			// We will run the 'chooser' function overSampleFraction more times than we project will be neccesary
			float overSampleFraction = 1.5f;

			foreach (var d in data)
			{
				var passRate = Math.Max(1e-9, (float)nPassed / (nTested + 1));
				var projectedGood = (total - nIterated) * passRate;
				var keepFrac = (target - nReturned) / projectedGood;
				var tryFrac = keepFrac*overSampleFraction;
				var innerFrac = keepFrac > 1.0 ? 1.0 : keepFrac/tryFrac;

				if (r.NextDouble() < tryFrac)
				{
					var p = chooser(d);
					nTested++;
					if (p != null)
					{
						nPassed++;
					}
					else
					{
						continue;
					}

					if (r.NextDouble() < innerFrac)
					{
						nReturned++;
						yield return p;
					}
				}

				nIterated++;
			}
            */
		}

	}
}
