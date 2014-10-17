using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.Linq;
using System.Timers;
using Bio.Algorithms.Assembly.Graph;
using Bio.Algorithms.Assembly.Padena.Scaffold;
using Bio.Util;
using Bio.Algorithms.Kmer;

namespace Bio.Algorithms.Assembly.Padena
{
    /// <summary>
    /// Implements a de bruijn based approach for
    /// assembly of DNA sequences.
    /// </summary>

    public class ParallelDeNovoAssembler : IDisposable
    {
        #region Fields
        /// <summary>
        /// Holds time interval between two progress events.
        /// </summary>
        private const int ProgressTimerInterval = 5 * 60 * 1000;

        /// <summary>
        /// User Input Parameter
        /// Length of k-mer.
        /// </summary>
        protected int _kmerLength;

        /// <summary>
        /// Timer to report progress.
        /// </summary>
        private Timer _progressTimer;

        /// <summary>
        /// Holds the current step number being executed.
        /// </summary>
        private int _currentStep;

        /// <summary>
        /// Flag to indicate whether graph building completed progress message is written to Console or not.
        /// </summary>
        private bool _graphBuildCompleted;

        /// <summary>
        /// Flag to indicate whether generate link completed progress message is written to Console or not.
        /// </summary>
        private bool _linkGenerationCompleted;

        #endregion

        /// <summary>
        /// Initializes a new instance of the ParallelDeNovoAssembler class.
        /// Sets thresholds to default values.
        /// Also initializes instances implementing different steps.
        /// </summary>
        public ParallelDeNovoAssembler()
        {
            ContigCoverageThreshold = -1;
            ErosionThreshold = -1;
            AllowErosion = false;
            // Initialize to default here.
            // Values set to -1 here will be reset based on input sequences.
            this._kmerLength = -1;
            this.DanglingLinksThreshold = -1;
            this.RedundantPathLengthThreshold = -1;

            // Contig and scaffold Builder are required modules. Set this to default.
            this.ContigBuilder = new SimplePathContigBuilder();

            // Default values for parameters used in building scaffolds.
            this.ScaffoldRedundancy = 2;
            this.Depth = 10;
            this.AllowKmerLengthEstimation = true;
        }

        #region Properties
        /// <summary>
        /// Provides the status to the subscribers.
        /// </summary>
        public static event EventHandler<StatusChangedEventArgs> StatusChanged;

        /// <summary>
        /// Gets the name of the current assembly algorithm used.
        /// This property returns the Name of our assembly algorithm i.e 
        /// Parallel De Novo algorithm.
        /// </summary>
        public string Name
        {
			get { return "Parallel de novo assembler"; }
        }

       

        /// <summary>
        /// Gets or sets the kmer length.
        /// </summary>
        public int KmerLength
        {
            get
            {
                return this._kmerLength;
            }

            set
            {
                if (value % 2 == 0)
                {
                    throw new ArgumentException("Cannot set the k-mer length to an even number.", "KmerLength");
                }                    
                this._kmerLength = value;
                this.AllowKmerLengthEstimation = false;
            }
        }

        /// <summary>
        /// Gets or sets a value indicating whether to estimate kmer length.
        /// </summary>
        public bool AllowKmerLengthEstimation { get; set; }

        /// <summary>
        /// Gets the assembler de-bruijn graph.
        /// </summary>
        public DeBruijnGraph Graph { get; protected set; }

        /// <summary>
        /// Gets or sets the instance that implements
        /// dangling links purging step.
        /// </summary>
        public IGraphErrorPurger DanglingLinksPurger { get; set; }

        /// <summary>
        /// Gets or sets the threshold length 
        /// for dangling link purger.
        /// </summary>
        [OutputAttribute]
        public int DanglingLinksThreshold { get; set; }

        /// <summary>
        /// Gets or sets a value indicating whether to allow erosion of the graph.
        /// </summary>
        [OutputAttribute]
        public bool AllowErosion { get; set; }

        /// <summary>
        /// Gets or sets the threshold length for eroding low coverage graph 
        /// ends. In case erosion step is not to be done, set this to 0.
        /// As an performance optimization in assembler process, erosion and 
        /// dangling link purging step are done together in a single step. 
        /// Note that because of this optimization, unless the danglingLinkPurger 
        /// implements IGraphErodePurger, erosion will not be done irrespective 
        /// of the threshold value provided. 
        /// </summary>
        [OutputAttribute]
        public int ErosionThreshold { get; set; }

        /// <summary>
        /// Gets or sets the instance that implements
        /// redundant paths purging step.
        /// </summary>
        public IGraphErrorPurger RedundantPathsPurger { get; set; }

        /// <summary>
        /// Gets or sets the length threshold 
        /// for redundant paths purger.
        /// </summary>
        [OutputAttribute]
        public int RedundantPathLengthThreshold { get; set; }

        /// <summary>
        /// Gets or sets instance of class implementing Low coverage contig removal.
        /// </summary>
        public ILowCoverageContigPurger LowCoverageContigPurger { get; set; }

        /// <summary>
        /// Gets or sets a value indicating whether to enable removal of low coverage contigs.
        /// </summary>
        [OutputAttribute]
        public bool AllowLowCoverageContigRemoval { get; set; }

        /// <summary>
        /// Gets or sets Threshold used for removing low-coverage contigs.
        /// </summary>
        [OutputAttribute]
        public double ContigCoverageThreshold { get; set; }

        /// <summary>
        /// Gets or sets the instance that implements
        /// contig building step.
        /// </summary>
        public IContigBuilder ContigBuilder { get; set; }

        /// <summary>
        /// Gets or sets the instance that implements
        /// scaffold building step.
        /// </summary>
        public IGraphScaffoldBuilder ScaffoldBuilder { get; set; }

        /// <summary>
        /// Gets or sets value of redundancy for building scaffolds.
        /// </summary>
        public int ScaffoldRedundancy { get; set; }

        /// <summary>
        /// Gets or sets the Depth for graph traversal in scaffold builder step.
        /// </summary>
        public int Depth { get; set; }

 
        #endregion

        /// <summary>
        /// For optimal graph formation, k-mer length should not be less 
        /// than half the length of the longest input sequence and 
        /// cannot be more than the length of the shortest input sequence. 
        /// Reference for estimating kmerlength from reads: Supplement material from 
        /// publication "ABySS: A parallel assembler for short read sequence data".
        /// </summary>
        /// <param name="sequences">List of input sequences.</param>
        /// <returns>Estimated optimal kmer length.</returns>
        public static int EstimateKmerLength(IEnumerable<ISequence> sequences)
        {
            if (sequences == null)
                throw new ArgumentNullException("sequences");

            //if (!Alphabets.CheckIsFromSameBase(sequences.First().Alphabet, Alphabets.DNA))
            //    throw new InvalidOperationException(Properties.Resource.CannotAssembleSequenceType);

            long minSeqLength = long.MaxValue, maxSeqLength = 0;

            // Get the min/max ranges for the sequences
            foreach (ISequence seq in sequences)
            {
                long seqCount = seq.Count;
                if (minSeqLength > seqCount)
                    minSeqLength = seqCount;
                if (maxSeqLength < seqCount)
                    maxSeqLength = seqCount;
            }

            // for optimal purpose, kmer length should be more than half of longest sequence
            float minLengthOfKmer = Math.Max(1, maxSeqLength / 2);
            float maxLengthOfKmer = minSeqLength;

            int kmerLength = minLengthOfKmer < maxLengthOfKmer
                                 ? (int) Math.Ceiling((minLengthOfKmer + maxLengthOfKmer)/2)
                                 : (int) Math.Floor(maxLengthOfKmer);

            // Make the kmer odd to avoid palindromes.
            if (kmerLength%2 == 0)
            {
                kmerLength++;
                if (kmerLength > maxLengthOfKmer)
                    kmerLength -= 2;
                if (kmerLength <= 0)
                    kmerLength = 1;
            }

            // Final sanity checks.
            if (maxLengthOfKmer < kmerLength)
				throw new InvalidOperationException("K-mer length was set too high");
            if (kmerLength <= 0)
				throw new InvalidOperationException("K-mer length must be greater than 0");

            // Bound to our max size based on data handling.
            return kmerLength > DeBruijnGraph.MaxKmerLength ? DeBruijnGraph.MaxKmerLength : kmerLength;
        }

        /// <summary>
        /// Assemble the list of sequence reads.
        /// </summary>
        /// <param name="inputSequences">List of input sequences.</param>
        /// <returns>Assembled output.</returns>
		public virtual PadenaAssembly Assemble(IEnumerable<ISequence> inputSequences)
        {
            if (inputSequences == null)
            {
                throw new ArgumentNullException("inputSequences");
            }

            // Remove ambiguous reads and set up fields for assembler process
            this.Initialize();

            // Step 1, 2: Create k-mers from reads and build de bruijn graph
            Stopwatch sw = Stopwatch.StartNew();
            this.CreateGraphStarted();
            this.CreateGraph(inputSequences);
            sw.Stop();

            this.CreateGraphEnded();
            this.TaskTimeSpanReport(sw.Elapsed);
            this.NodeCountReport();
           
            // Estimate and set default value for erosion and coverage thresholds
            sw = Stopwatch.StartNew();
            this.EstimateDefaultValuesStarted();
            this.EstimateDefaultThresholds();
            sw.Stop();

            this.EstimateDefaultValuesEnded();
            this.TaskTimeSpanReport(sw.Elapsed);
            
            // Step 3: Remove dangling links from graph
            sw = Stopwatch.StartNew();
            this.UndangleGraphStarted();
            this.UnDangleGraph();
            sw.Stop();

            this.UndangleGraphEnded();
            this.TaskTimeSpanReport(sw.Elapsed);
            this.NodeCountReport();
            
            // Step 4: Remove redundant paths from graph
            sw = Stopwatch.StartNew();
            this.RemoveRedundancyStarted();
            this.RemoveRedundancy();
            this.NodeCountReport();
            
            // Perform dangling link purger step once more.
            // This is done to remove any links created by redundant paths purger.
			RaiseMessage(string.Format(CultureInfo.CurrentCulture, "Second round undangle graph start.", DateTime.Now));
            this.UnDangleGraph();
			RaiseMessage(string.Format(CultureInfo.CurrentCulture, "Second round undangle graph end.", DateTime.Now));
            
            // Report end after undangle
            sw.Stop();
            this.RemoveRedundancyEnded();
            this.TaskTimeSpanReport(sw.Elapsed);
            this.NodeCountReport();

            // Step 5: Build Contigs
            sw = Stopwatch.StartNew();
            this.BuildContigsStarted();
            IEnumerable<ISequence> contigSequences = this.BuildContigs();
            sw.Stop();

            this.BuildContigsEnded();
            this.TaskTimeSpanReport(sw.Elapsed);
            
            PadenaAssembly result = new PadenaAssembly();
            result.AddContigs(contigSequences);

            return result;
        }       

        /// <summary>
        /// Implements dispose to suppress GC finalize
        /// This is done as one of the methods uses ReadWriterLockSlim
        /// which extends IDisposable.
        /// </summary>
        public void Dispose()
        {
            this.Dispose(true);
            GC.SuppressFinalize(this);
        }

        /// <summary>
        /// Estimates and sets erosion and coverage threshold for contigs.
        /// Median value of kmer coverage is set as default value.
        /// Reference: ABySS Release Notes 1.1.1 - "The default threshold 
        /// is the square root of the median k-mer coverage".
        /// </summary>
        protected virtual void EstimateDefaultThresholds()
        {
            if (this.AllowErosion || this.AllowLowCoverageContigRemoval)
            {
                // In case of low coverage data, set default as 2.
                // Reference: ABySS Release Notes 1.0.15
                // Before calculating median, discard thresholds less than 2.
                List<long> kmerCoverage = this.Graph.GetNodes().AsParallel().Aggregate(
                    new List<long>(),
                    (kmerList, n) =>
                    {
                        if (n.KmerCount > 2)
                        {
                            kmerList.Add(n.KmerCount);
                        }

                        return kmerList;
                    });

                double threshold;
                if (kmerCoverage.Count == 0)
                {
                    threshold = 2; // For low coverage data, set default as 2
                }
                else
                {
                    kmerCoverage.Sort();
                    int midPoint = kmerCoverage.Count / 2;
                    double median = (kmerCoverage.Count % 2 == 1 || midPoint == 0) ?
                        kmerCoverage[midPoint] :
                        ((float)(kmerCoverage[midPoint] + kmerCoverage[midPoint - 1])) / 2;
                    
                    threshold = Math.Sqrt(median);
                }

                // Set coverage threshold
                if (this.AllowLowCoverageContigRemoval && this.ContigCoverageThreshold == -1)
                {
                    this.ContigCoverageThreshold = threshold;
                }

                if (this.AllowErosion && this.ErosionThreshold == -1)
                {
                    // Erosion threshold is an int, so round it off
                    this.ErosionThreshold = (int)Math.Round(threshold);
                }
            }
        }

        /// <summary>
        /// Step 1: Building k-mers from sequence reads
        /// Step 2: Build de bruijn graph for input set of k-mers.
        /// Sets the _assemblerGraph field.
        /// </summary>
        protected virtual void CreateGraph(IEnumerable<ISequence> inputSequences)
        {
            this.Graph = new DeBruijnGraph(this._kmerLength);
            this.Graph.Build(inputSequences);

            // Recapture the kmer length to keep them in sync.
            this._kmerLength = this.Graph.KmerLength;
        }

        /// <summary>
        /// Step 3: Remove dangling links from graph.
        /// </summary>
        protected void UnDangleGraph()
        {
            if (this.DanglingLinksPurger != null && this.DanglingLinksThreshold > 0)
            {
                DeBruijnPathList danglingNodes = null;

                // Observe lengths of dangling links in the graph
                // This is an optimization - instead of incrementing threshold by 1 and 
                // running the purger iteratively, we first determine the lengths of the 
                // danglings links found in the graph and run purger only for those lengths.
                this.DanglingLinksPurger.LengthThreshold = this.DanglingLinksThreshold - 1;

                IEnumerable<int> danglingLengths;
                IGraphEndsEroder graphEndsEroder = this.DanglingLinksPurger as IGraphEndsEroder;
                if (graphEndsEroder != null && this.AllowErosion)
                {
                    // If eroder is implemented, while getting lengths of dangling links, 
                    // it also erodes the low coverage ends, this marks any node for deletion below a threshold.

                    //TODO: Verify that this does enumerate all dangling ends, the concern is that if a dangling end of length 7 and 2
                    //arrive at a node which itself would be of dangling node of length 2 without these "dangling ends" then a dangling end of 9
                    // (which it would be without either the 7 or 2 end) might not be reported.
                    danglingLengths = graphEndsEroder.ErodeGraphEnds(this.Graph, this.ErosionThreshold);
                }
                else
                {
                    // Perform dangling purger at all incremental values till dangleThreshold.
                    danglingLengths = Enumerable.Range(1, this.DanglingLinksThreshold - 1);
                }

                // Erosion is to be only once. Reset erode threshold to -1.
                this.ErosionThreshold = -1;
              
            
                // Start removing dangling links
                foreach (int threshold in danglingLengths)
                {
                    if (this.Graph.NodeCount >= threshold)
                    {
                        this.DanglingLinksPurger.LengthThreshold = threshold;
                        danglingNodes = this.DanglingLinksPurger.DetectErroneousNodes(this.Graph);
                        this.DanglingLinksPurger.RemoveErroneousNodes(this.Graph, danglingNodes);
                    }
                }

                // Removing dangling links can in turn create more dangling links
                // In order to remove all links within threshold, we therefore run
                // purger at threshold length until there is no more change in graph.
                do
                {
                    danglingNodes = null;
                    if (this.Graph.NodeCount >= this.DanglingLinksThreshold)
                    {
                        this.DanglingLinksPurger.LengthThreshold = this.DanglingLinksThreshold;
                        danglingNodes = this.DanglingLinksPurger.DetectErroneousNodes(this.Graph);
                        this.DanglingLinksPurger.RemoveErroneousNodes(this.Graph, danglingNodes);
                    }
                }
                while (danglingNodes != null && danglingNodes.Paths.Count > 0);
            }
        }

        /// <summary>
        /// Step 4: Remove redundant paths from graph.
        /// </summary>
        protected void RemoveRedundancy()
        {
            if (this.RedundantPathsPurger != null)
            {
                DeBruijnPathList redundantNodes;
                do
                {
                    redundantNodes = this.RedundantPathsPurger.DetectErroneousNodes(this.Graph);
                    this.RedundantPathsPurger.RemoveErroneousNodes(this.Graph, redundantNodes);
                }
                while (redundantNodes.Paths.Count > 0);
            }
        }

        /// <summary>
        /// Step 5: Build contigs from de bruijn graph.
        /// If coverage threshold is set, remove low coverage contigs.
        /// </summary>
        /// <returns>List of contig sequences.</returns>
        protected virtual IEnumerable<ISequence> BuildContigs()
        {
            if (this.ContigBuilder == null)
            {
				throw new InvalidOperationException("nullcontigbuilder");
            }

            // Step 5.1: Remove low coverage contigs
            if (this.AllowLowCoverageContigRemoval && this.ContigCoverageThreshold > 0)
            {
                this.LowCoverageContigPurger.RemoveLowCoverageContigs(this.Graph, this.ContigCoverageThreshold);
            }

            // Step 5.2: Build Contigs
            return this.ContigBuilder.Build(this.Graph);
        }

        /// <summary>
        /// Dispose field instances.
        /// </summary>
        /// <param name="disposeManaged">If disposeManaged equals true, clean all resources.</param>
        protected virtual void Dispose(bool disposeManaged)
        {
            if (disposeManaged)
            {
                if (this.ScaffoldBuilder != null)
                {
                    this.ScaffoldBuilder.Dispose();
                }

                this.Graph = null;
                this.DanglingLinksPurger = null;
                this.RedundantPathsPurger = null;
                this.ContigBuilder = null;
                this.ScaffoldBuilder = null;

                if (this._progressTimer != null)
                {
                    this._progressTimer.Dispose();
                    this._progressTimer = null;
                }
            }
        }

        /// <summary>
        /// Sets up fields for the assembly process.
        /// </summary>
        protected void Initialize()
        {
            this._currentStep = 0;
            this._progressTimer = new Timer(ProgressTimerInterval);
            this._progressTimer.Elapsed += this.ProgressTimerElapsed;

            //this.RaiseMessage((CultureInfo.CurrentCulture, "Initializing - Start time: {0:yyyy-MM-dd-HH:mm:ss.fff}", DateTime.Now));

            // Reset parameters not set by user, based on sequenceReads
            if (this._kmerLength <= 0)
            {
                throw new InvalidOperationException("Kmer length must be less than the length of the shortest sequence, and should be greater than half the length of the longest sequence.");
            }

            // TODO: Force an odd kmer length to avoid palindromes
            // Need to evaluate this more - for now rely on the user to
            // not pass in bad data.
            if (_kmerLength % 2 == 0)
                _kmerLength++;

            // Enforce our boundaries (same as DeBruijnGraph code)
            _kmerLength = Math.Max(1, Math.Min(DeBruijnGraph.MaxKmerLength, _kmerLength));

            if (this.DanglingLinksThreshold == -1)
            {
                this.DanglingLinksThreshold = this._kmerLength + 1;
            }

            if (this.RedundantPathLengthThreshold == -1)
            {
                // Reference for default threshold for redundant path purger:
                // ABySS Release Notes 1.1.2 - "Pop bubbles shorter than N bp. The default is b=3*(k + 1)."
                this.RedundantPathLengthThreshold = 3 * (this._kmerLength + 1);
            }

            this.InitializeDefaultGraphModifiers();

            //RaiseMessage(string.Format(CultureInfo.CurrentCulture, "Initializing - End time: {0:yyyy-MM-dd-HH:mm:ss.fff}", DateTime.Now));
        }

        /// <summary>
        /// Initializes the above defined fields. For each step in assembly
        /// we use a separate class for implementation. This method assigns 
        /// these variables to classes with desired implementation.
        /// </summary>
        protected void InitializeDefaultGraphModifiers()
        {
            // Assign uninitialized fields to default values
            if (this.DanglingLinksPurger == null)
            {
                this.DanglingLinksPurger = new DanglingLinksPurger();
            }

            if (this.RedundantPathsPurger == null)
            {
                this.RedundantPathsPurger = new RedundantPathsPurger(this.RedundantPathLengthThreshold);
            }

            if (this.LowCoverageContigPurger == null)
            {
                this.LowCoverageContigPurger = new SimplePathContigBuilder();
            }
        }


        #region ProgressReportMethodsAndProperties

        /// <summary>
        /// Raises status event.
        /// </summary>
        [System.Diagnostics.CodeAnalysis.SuppressMessage("Microsoft.Design", "CA1030:UseEventsWhereAppropriate")]
        public static void RaiseStatusEvent(string msg)
        {
			//Console.WriteLine(_statusMessage);
            if (StatusChanged != null)
            {
                StatusChanged.Invoke(null, new StatusChangedEventArgs(msg));
            }
        }

        /// <summary>
        /// Method to handle ProgressTimer elapsed event.
        /// </summary>
        /// <param name="sender">Progress timer.</param>
        /// <param name="e">Event arguments.</param>
        protected void ProgressTimerElapsed(object sender, ElapsedEventArgs e)
        {
            switch (this._currentStep)
            {
                case 2:
                    if (!this.Graph.GraphBuildCompleted)
                    {
                        if (this.Graph.SkippedSequencesCount > 0)
                        {
						 RaiseMessage(string.Format(CultureInfo.CurrentCulture, "Creating Graph Progress: {0} sequence(s) skipped out of {1} processed sequences.  Currently {2} Nodes in Graph", this.Graph.SkippedSequencesCount, this.Graph.ProcessedSequencesCount,this.Graph.NodeCount)); 
                        }
                        else
                        {
						    RaiseMessage(string.Format(CultureInfo.CurrentCulture, "{0} sequence(s) processed.", this.Graph.ProcessedSequencesCount));
                        }
                    }
                    else
                    {
                        if (!this._graphBuildCompleted)
                        {
                            this._graphBuildCompleted = this.Graph.GraphBuildCompleted;
						    RaiseMessage(string.Format(CultureInfo.CurrentCulture,"   Graph built successfully - Processed {0} sequences.", this.Graph.ProcessedSequencesCount));
						    RaiseMessage(string.Format(CultureInfo.CurrentCulture,"   Generate Links Started."));
                        }
                        else
                        {
                            if (!this._linkGenerationCompleted && this.Graph.LinkGenerationCompleted)
                            {
                                this._linkGenerationCompleted = this.Graph.LinkGenerationCompleted;
                                RaiseMessage(string.Format(CultureInfo.CurrentCulture,
								"   Generate Links Ended."));
                            }
                            else
                            {
							     RaiseMessage(".");
                            }
                        }
                    }
                    break;
                default:
				    RaiseMessage( ".");
                    break;
            }
            var mem=Bio.Util.Logging.OutputInformation.GetMemoryUsage();
            RaiseMessage("Currently using: " + mem.PeakWorkingSet + " working set memory.)");
            //if (mem.NumericWorkingSet > 5e9)
           // {
           //    RaiseMessage("Manually invoking GC");
           //     GC.Collect();
           //     mem=Bio.Util.Logging.OutputInformation.GetMemoryUsage();
           //     RaiseMessage("Currently using: " + mem.PeakWorkingSet + " working set memory.");
           // }
            
        }

        /// <summary>
        /// Raises status changed event with Graph creating started status message.
        /// </summary>
        protected void CreateGraphStarted()
        {
			RaiseMessage(string.Format(CultureInfo.CurrentCulture, "Step 1 & 2: Create Kmer and Graph - Start time: {0:yyyy-MM-dd-HH:mm:ss.fff}", DateTime.Now));
            this._currentStep = 2;
            this._progressTimer.Start();
        }

        /// <summary>
        /// Raises status changed event with Graph creating ended status message.
        /// </summary>
        protected void CreateGraphEnded()
        {
            this._progressTimer.Stop();

            if (!this._graphBuildCompleted)
            {
                this._graphBuildCompleted = this.Graph.GraphBuildCompleted;
				RaiseMessage(string.Format(CultureInfo.CurrentCulture,"   Graph built successfully - Processed {0} sequences.", this.Graph.ProcessedSequencesCount));
                RaiseMessage(string.Format(CultureInfo.CurrentCulture, "   Generate Links Started."));
            }

            if (!this._linkGenerationCompleted && this.Graph.LinkGenerationCompleted)
            {
                this._linkGenerationCompleted = this.Graph.LinkGenerationCompleted;
				RaiseMessage(string.Format(CultureInfo.CurrentCulture, "   Generate Links Ended."));
            }

			RaiseMessage(string.Format(CultureInfo.CurrentCulture, "Step 1 & 2: Create Kmer and Graph - End time: {0:yyyy-MM-dd-HH:mm:ss.fff}", DateTime.Now));
        }

        /// <summary>
        /// Report the time a task took
        /// </summary>
        /// <param name="ts"></param>
        protected void TaskTimeSpanReport(TimeSpan ts)
        {
			RaiseMessage(string.Format(CultureInfo.CurrentCulture, "Total Time for Task: {0:c}", ts));            
        }

        /// <summary>
        /// Raise event to report the number of nodes currently in graph.
        /// </summary>
        protected void NodeCountReport()
        {
            RaiseMessage(string.Format(CultureInfo.CurrentCulture, "Nodes in graph: {0}", this.Graph.NodeCount));
        }

        /// <summary>
        /// Raises status changed event with EstimateDefaultValues started status message.
        /// </summary>
        protected void EstimateDefaultValuesStarted()
        {
			RaiseMessage(string.Format(CultureInfo.CurrentCulture, "Estimating default values - Start time: {0:yyyy-MM-dd-HH:mm:ss.fff}", DateTime.Now));
        }

        /// <summary>
        /// Raises status changed event with EstimateDefaultValues ended status message.
        /// </summary>
        protected void EstimateDefaultValuesEnded()
        {
			RaiseMessage(string.Format(CultureInfo.CurrentCulture, "Estimating default values - End time: {0:yyyy-MM-dd-HH:mm:ss.fff}", DateTime.Now));
        }

        /// <summary>
        /// Raises status changed event with UndangleGraph started status message.
        /// </summary>
        protected void UndangleGraphStarted()
        {
            RaiseMessage(string.Format(CultureInfo.CurrentCulture, "UndangleGraph - Start time: {0:yyyy-MM-dd-HH:mm:ss.fff}", DateTime.Now));
            this._currentStep = 3;
            this._progressTimer.Start();
        }

        /// <summary>
        /// Raises status changed event with UndangleGraph ended status message.
        /// </summary>
        protected void UndangleGraphEnded()
        {
            this._progressTimer.Stop();
			RaiseMessage(string.Format(CultureInfo.CurrentCulture, "UndangleGraph - End time: {0:yyyy-MM-dd-HH:mm:ss.fff}", DateTime.Now));
        }

        /// <summary>
        /// Raises status changed event with RemoveRedundancy started status message.
        /// </summary>
        protected void RemoveRedundancyStarted()
        {
            RaiseMessage(string.Format(CultureInfo.CurrentCulture, "RemoveRedundancy - Start time: {0:yyyy-MM-dd-HH:mm:ss.fff}", DateTime.Now));
            RaiseMessage("Redundant path threshold is: " + this.RedundantPathLengthThreshold.ToString());
        }

        /// <summary>
        /// Raises status changed event with RemoveRedundancy ended status message.
        /// </summary>
        protected void RemoveRedundancyEnded()
        {
            this._progressTimer.Stop();
			RaiseMessage(string.Format(CultureInfo.CurrentCulture, "End Redundancy Removal Step", DateTime.Now));
        }

        /// <summary>
        /// Raises status changed event with BuildContigs started status message.
        /// </summary>
        protected void BuildContigsStarted()
        {
			RaiseMessage(string.Format(CultureInfo.CurrentCulture, "Start building Contigs", DateTime.Now));
            this._currentStep = 5;
            this._progressTimer.Start();
        }

        /// <summary>
        /// Raises status changed event with BuildContigs ended status message.
        /// </summary>
        protected void BuildContigsEnded()
        {
            this._progressTimer.Stop();
			RaiseMessage(string.Format(CultureInfo.CurrentCulture, "End Building Contigs", DateTime.Now));
            this._currentStep = 0;
        }

        /// <summary>
        /// Raises status changed event with BuildScaffolds started status message.
        /// </summary>
        protected void BuildScaffoldsStarted()
        {
			RaiseMessage(string.Format(CultureInfo.CurrentCulture, "Start building scaffolds", DateTime.Now),true);
            this._currentStep = 6;
            this._progressTimer.Start();
        }

        /// <summary>
        /// Raises status changed event with BuildScaffolds ended status message.
        /// </summary>
        protected void BuildScaffoldsEnded()
        {
            this._progressTimer.Stop();
            this._currentStep = 0;
			RaiseMessage(string.Format(CultureInfo.CurrentCulture, "BuildScaffolds - End time: {0:yyyy-MM-dd-HH:mm:ss.fff}", DateTime.Now),true);
        }
        /// <summary>
        /// Raise a message
        /// </summary>
        /// <param name="msg"></param>
        public void RaiseMessage(string msg, bool indent = true)
        {
            if (indent)
            {
                msg = "\t" + msg;
            }
            RaiseStatusEvent(msg);
        }
        #endregion
    }
}
