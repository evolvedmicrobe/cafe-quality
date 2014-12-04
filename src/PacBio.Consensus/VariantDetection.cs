using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using PacBio.Align;
using PacBio.IO;
using PacBio.Utils;

namespace PacBio.Consensus
{
    public struct VariantDetectionConfig
    {
        public bool CountAmbiguous;
        public int WindowRadius;
        public double MinLocalAccuracy;
        public int MinCoverage;
        public int MaxIterations;
        public double BinomialTestThreshold;
        public double ConvergenceCriterion;
        public double ConfidenceInterval;
        public double MaxPValue;
        public double BetaPrior;
        public double MinErrorRate;
    }


    // LoFreqNQ model from:
    //
    // Wilm A, Aw PP, Bertrand D, Yeo GH, Ong SH, Wong CH, Khor CC, Petric R, Hibberd
    // ML, Nagarajan N. "LoFreq: a sequence-quality aware, ultra-sensitive variant caller
    // for uncovering cell-population heterogeneity from high-throughput sequencing
    // datasets." Nucleic Acids Res. 2012 Dec;40(22):11189-201. doi: 10.1093/nar/gks918. 
    // Epub 2012 Oct 12. PubMed PMID: 23066108; PubMed Central PMCID: PMC3526318.
    public class VariantDetection
    {
        private readonly VariantDetectionConfig _config;
        private double[,] _parameters;


        public VariantDetection(VariantDetectionConfig config)
        {
            _config = config;
            _parameters = null;
        }


       

        public struct SiteCount
        {
            public double[] Counts;
            public int TemplatePosition;
            public double Coverage;
            public int RefBase;
        }


        private enum Nucleotide
        {
            N = -1,
            A = 0,
            C = 1,
            G = 2,
            T = 3,
            Length = 4
        };


        private static Nucleotide NucToIndex(char nuc)
        {
            switch (nuc)
            {
                case 'A':
                    return Nucleotide.A;
                case 'C':
                    return Nucleotide.C;
                case 'G':
                    return Nucleotide.G;
                case 'T':
                    return Nucleotide.T;
            }
            return Nucleotide.N;
        }

        private static char IndexToNuc(Nucleotide nuc)
        {
            switch (nuc)
            {
                case Nucleotide.A:
                    return 'A';
                case Nucleotide.C:
                    return 'C';
                case Nucleotide.G:
                    return 'G';
                case Nucleotide.T:
                    return 'T';
            }
            return 'N';
        }


        private class Counter<T>
        {
            private Dictionary<T, int> dict;

            public Counter()
            {
                dict = new Dictionary<T, int>();
            }

            public Counter(IEnumerable<T> e)
            {
                dict = new Dictionary<T, int>();
                AddRange(e);
            }

            public void AddRange(IEnumerable<T> e)
            {
                foreach (var i in e)
                {
                    Add(i);
                }
            }

            public void Add(T i)
            {
                if (dict.ContainsKey(i))
                    dict[i] += 1;
                else
                    dict[i] = 1;
            }

            public IEnumerable<T> Keys { get { return dict.Keys; } }

            public IEnumerator<KeyValuePair<T, int>> GetEnumerator()
            {
                return dict.GetEnumerator();
            }
        }

        public static IEnumerable<Tuple<int, char>> ReadByTemplateIndex(IAlnSummary alnSummary)
        {
            var aln = alnSummary.ReadAlignmentStrings(Orientation.Genomic);
            var tplPos = alnSummary.TemplateStart;

            for (var i = 0; i < aln.Item1.Length; i++)
            {
                var refBase = aln.Item1[i];
                var qryBase = aln.Item2[i];

                if (qryBase != '-')
                    yield return Tuple.Create(tplPos, qryBase);

                if (refBase != '-')
                    tplPos += 1;
            }
        }


        public static bool SpansWindow(IAlnSummary alnSummary, int templateStart, int templateEnd)
        {
            return alnSummary.TemplateStart <= templateStart &&
                   alnSummary.TemplateEnd >= templateEnd;
        }


        public static string ClipRead(IEnumerable<Tuple<int, char>> seqByPos,
                                       int templateStart, int templateEnd)
        {
            var bases = seqByPos
                .Where(sp => sp.Item1 >= templateStart &&
                             sp.Item1 < templateEnd)
                .Select(sp => sp.Item2)
                .ToArray();
            return new string(bases);
        }


        private static SiteCount[] CountNucleotides(IEnumerable<IAlnSummary> alignments,
                                                    TrialTemplate template, int templateStart,
                                                    bool countAmbiguous, double minAccuracy,
                                                    int windowRadius)
        {
            const int nBases = (int) Nucleotide.Length;
            var siteCounts = new SiteCount[template.Length];
            var tplString = template.GetSequence(Strand.Forward).ToUpper();

            // initialize counts
            for (var s = 0; s < siteCounts.Length; s++)
            {
                siteCounts[s].Counts = new double[nBases];
                siteCounts[s].TemplatePosition = s;
            }

            // per-site locks, should be fine-grained enough
            var siteWindows = siteCounts.Length.Fill(_ => new Counter<string>());
            var siteLocks = siteCounts.Length.Fill(_ => new object());

            alignments.Select((alnSummary, idx) => Tuple.Create(alnSummary, idx)).ParForEach(
                alnSummaryIdx =>
                {
                    var alnSummary = alnSummaryIdx.Item1;
                    var idx = alnSummaryIdx.Item2;

                    var seqByPos = ReadByTemplateIndex(alnSummary).ToList();

            
                    for (int siteIdx = 0; siteIdx < template.Length; siteIdx++)
                    {
                        var windowStart = Math.Max(0, siteIdx - windowRadius) + templateStart;
                        var windowEnd = Math.Min(template.Length, siteIdx + windowRadius + 1) + templateStart;

                        if (!SpansWindow(alnSummary, windowStart, windowEnd))
                            continue;

                        var qry = ClipRead(seqByPos, windowStart, windowEnd).ToUpper();

                        lock (siteLocks[siteIdx])
                        {
                            siteWindows[siteIdx].Add(qry);
                        }
                    }
                });

            Enumerable.Range(0, siteCounts.Length).ParForEach(
                siteIdx =>
                {
                    var windowStart = Math.Max(0, siteIdx - windowRadius);
                    var windowEnd = Math.Min(template.Length, siteIdx + windowRadius + 1);
                    var tpl = tplString.Substring(windowStart, windowEnd - windowStart);

                    foreach (var seqCount in siteWindows[siteIdx])
                    {
                        var obs = LocalAlignment.CountRow(tpl, seqCount.Key, siteIdx - windowStart, minAccuracy);

                        if (!obs.Any())
                            continue;

                        if (obs.Count > 1 && !countAmbiguous)
                        {
                            //Log(LogLevel.INFO, "Skipping {0} alignments at site {1} due to ambiguity: {2}",
                            //    alnCount.Value, siteIdx + 1, String.Join(", ", obs));
                            continue;
                        }

                        lock (siteLocks[siteIdx])
                        {
                            foreach (var b in obs)
                            {
                                // TODO (lhepler) -- defensive programming here? LocalAlignment will only ever provide ACGT
                                var bIdx = (int) NucToIndex(b);
                                siteCounts[siteIdx].Counts[bIdx] += ((double) seqCount.Value)/obs.Count;
                            }
                        }
                    }
                });

            for (var s = 0; s < siteCounts.Length; s++)
            {
                siteCounts[s].Coverage = siteCounts[s].Counts.Sum();
                siteCounts[s].RefBase  = siteCounts[s].Counts.IMax();
            }

            return siteCounts;
        }


        private static double SiteAltBasePValue(SiteCount siteCount, int refBase, int altBase,
                                                double[,] parameters, double minRate)
        {
            const int nBases = (int) Nucleotide.Length;

            var n = siteCount.Coverage;
            var k = siteCount.Counts[altBase];

            // if we didn't observe anything, then we have nothing to report
            if (k <= 0)
                return 1.0;

            if (n <= k)
                return 0.0;

            if (refBase == altBase)
                return 0.0;

            /*
            // TODO: this should probably be a beta-binomial model,
            // but computing the 3F2 hypergeometric series is (ahem) nontrivial.
            // So approximate the dispersion by overestimating the rate
            // for the number of samples at this site
            var ciHigh  = 1 - threshold / (2 * n);  // TODO: appropriate Bonferroni correction?
            var alpha   = n * parameters[refBase, altBase];
            var beta    = n - alpha;
            var mean    = alpha / (alpha + beta);
            var rate    = MathUtils.FindRootApprox(x => MathUtils.BetaRegularized(alpha, beta, x) - ciHigh,
                                                   mean, 1);
             */
            var rate = Enumerable.Range(0, nBases).Select(b => parameters[refBase, b]).Sum();
            var invRate = 1.0 - Math.Max(minRate, rate);

            // Binomial Survival Function (inclusive)
            return 1.0 - MathUtils.BetaRegularized(n - k + 1, k, invRate);
        }


        private static void Expect(IList<SiteCount> siteCounts, double[,] parameters,
                                   bool[,] isError, double threshold, double minRate)
        {
            const int nBases = (int) Nucleotide.Length;
            var nSites = siteCounts.Count;

            isError.Fill(true);

            for (var site = 0; site < nSites; site++)
            {
                var refBase = siteCounts[site].RefBase;

                isError[site, refBase] = false;

                for (var altBase = 0; altBase < nBases; altBase++)
                {
                    if (altBase == refBase)
                        continue;

                    var pValue = SiteAltBasePValue(siteCounts[site], refBase, altBase,
                                                   parameters, minRate);

                    // TODO: Bonferroni correction? (threshold/n?)
                    if (pValue < threshold)
                    {
                        isError[site, altBase] = false;
                    }
                }
            }
        }


        private static void Maximize(IList<SiteCount> siteCounts, bool[,] isError,
                                     double[,] parameters, double[,] tmp, int minCoverage)
        {
            const int nBases = (int) Nucleotide.Length;
            var nSites = siteCounts.Count;

            parameters.Fill(0);
            tmp.Fill(0);

            // transition probability from X to Y P(X -> Y)
            // is maximized by the following, where C is the column, Z_C \in {R,V} is the model
            // of the column, where model R has only a reference base and 12 parameters
            // for sequencing error and model V is variant, b(C) is the reference base at column C,
            // and n_X^C is the count of nucleotide X at column C:
            // P(X -> Y) = (sum_{Z_C=R,b(C)=X} n_Y^C) / (sum_{Z_C=R,b(C)=X} n_X^C)
            for (var site = 0; site < nSites; site++)
            {
                var isVariant = Enumerable.Range(0, nBases)
                    .Aggregate(0, (acc, b) => isError[site, b] ? 0 : 1) > 1;

                if (isVariant || siteCounts[site].Coverage < minCoverage)
                    continue;

                var refBase = siteCounts[site].RefBase;

                for (var altBase = 0; altBase < nBases; altBase++)
                {
                    tmp[refBase, altBase] += siteCounts[site].Counts[altBase];
                }
            }

            for (var from = 0; from < nBases; from++)
            {
                for (var to = 0; to < nBases; to++)
                {
                    if (from == to)
                        continue;

                    if (tmp[from, from] > 0)
                    {
                        // if we can divide, do so
                        parameters[from, to] = tmp[from, to]/tmp[from, from];
                    }
                    else if (tmp[from, to] > 0)
                    {
                        // otherwise if we're 0/0, then 0
                        parameters[from, to] = 1.0;
                    }
                    else
                    {
                        // otherwise we're completely in error
                        parameters[from, to] = 0.0;
                    }
                }
            }
        }


        private static bool IsConverged(Tuple<double[,], double[,]> mats, double criterion)
        {
            var dim1 = mats.Item1.GetLength(0);
            var dim2 = mats.Item1.GetLength(1);

            Debug.Assert(dim1 == mats.Item2.GetLength(0));
            Debug.Assert(dim2 == mats.Item2.GetLength(1));

            return Enumerable.Range(0, dim1)
                .SelectMany(
                    i =>
                        Enumerable.Range(0, dim2)
                        .Where(j => i != j)
                        .Select(j => Math.Abs(mats.Item1[i, j] - mats.Item2[i, j]) < criterion))
                .Aggregate(true, (x, y) => x && y);
        }


        private static IEnumerable<Variant> ComputeVariants(IList<SiteCount> siteCounts,
                                                            double[,] parameters, bool[,] isError,
                                                            TrialTemplate template, int templateStart,
                                                            VariantDetectionConfig config)
        {
            const int nBases = (int) Nucleotide.Length;

            // confidence intervals
            var ciLow  = (1 - config.ConfidenceInterval) / 2;
            var ciHigh = 1 - ciLow;

            var tplString = template.GetSequence(Strand.Forward);

            for (var site = 0; site < siteCounts.Count; site++)
            {
                var refBase = siteCounts[site].RefBase;

                for (var altBase = 0; altBase < nBases; altBase++)
                {
                    var tplPos = siteCounts[site].TemplatePosition;
                    var tplBaseChar = tplString[tplPos];
                    var altBaseChar = IndexToNuc((Nucleotide) altBase);
                    var pValue = SiteAltBasePValue(siteCounts[site], refBase, altBase, parameters, config.MinErrorRate);

                    if (pValue > config.MaxPValue ||
                        siteCounts[site].Coverage < config.MinCoverage ||
                        altBaseChar == tplBaseChar)
                        continue;

                    var k = siteCounts[site].Counts[altBase];
                    var n = Enumerable.Range(0, nBases)
                        .Where(b => !isError[site, b])
                        .Select(b => siteCounts[site].Counts[b])
                        .Sum();

                    var alpha = config.BetaPrior + k;
                    var beta  = config.BetaPrior + n - k;

                    if (alpha < 0 || beta < 0)
                        throw new ApplicationException(String.Format("Invalid alpha ({0}) and beta ({1})",
                                                                     alpha, beta));

                    var mean = alpha/(alpha + beta);

                    // if n == k then 100% have this mutation
                    var fractionLow = k <= 0 ? 0.0 : n <= k ? 1.0
                        : MathUtils.FindRootApprox(x => MathUtils.BetaRegularized(alpha, beta, x) - ciLow,
                                                   0, mean);

                    var fractionHigh = k <= 0 ? 0.0 : n <= k ? 1.0
                        : MathUtils.FindRootApprox(x => MathUtils.BetaRegularized(alpha, beta, x) - ciHigh,
                                                   mean, 1);

                    yield return new Variant
                                 {
                                     TemplatePosition = tplPos + templateStart,
                                     RefBase = tplBaseChar,
                                     AltBase = altBaseChar,
                                     Fraction = mean,
                                     FractionLow = fractionLow,
                                     FractionHigh = fractionHigh,
                                     PValue = pValue,
                                     PseudoCounts = k,
                                     Coverage = siteCounts[site].Coverage
                                 };
                }
            }
        }


        public IEnumerable<Variant> LearnApply(
            IList<IAlnSummary> alignments,
            TrialTemplate template, int templateStart)
        {
            return Run(alignments, template, templateStart, _config, out _parameters);
        }


        public bool IsLearned()
        {
            return _parameters != null;
        }


        public IEnumerable<Variant> Apply(
            IList<IAlnSummary> alignments,
            TrialTemplate template, int templateStart)
        {
            const int nBases = (int) Nucleotide.Length;

            var siteCounts = CountNucleotides(alignments,
                                              template, templateStart,
                                              _config.CountAmbiguous, _config.MinLocalAccuracy,
                                              _config.WindowRadius);
            var nSites = siteCounts.Length;
            var isError = new bool[nSites, nBases];
            
            Expect(siteCounts, _parameters, isError, _config.BinomialTestThreshold, _config.MinErrorRate);

            return ComputeVariants(siteCounts, _parameters, isError, template, templateStart, _config);
        }


        public static IEnumerable<Variant> Run(IList<IAlnSummary> alignments,
                                               TrialTemplate template, int templateStart,
                                               VariantDetectionConfig config,
                                               out double[,] finalParameters)
        {
            const int nBases = (int) Nucleotide.Length;

            var siteCounts = CountNucleotides(alignments,
                                              template, templateStart,
                                              config.CountAmbiguous, config.MinLocalAccuracy,
                                              config.WindowRadius);

            var nSites = siteCounts.Length;
            var isError = new bool[nSites, nBases];
            var parameters = new Tuple<double[,], double[,]> (new double[nBases, nBases],
                                                              new double[nBases, nBases]);
            var tmp = new double[nBases, nBases];

            // initialize isError (everything is an error but plurality base at each site)
            isError.Fill((s, b) => b != siteCounts[s].RefBase);

            // initialize the parameters 
            Maximize(siteCounts, isError, parameters.Item1, tmp, config.MinCoverage);

            var converged = false;

            // this needs a default value
            finalParameters = null;

            for (var iter = 0; iter < config.MaxIterations; iter++)
            {
                // alternate between parameter sets, start with src = Item1
                var srcParameters = iter%2 == 0 ? parameters.Item1 : parameters.Item2;
                var dstParameters = iter%2 == 0 ? parameters.Item2 : parameters.Item1;

           
                // Except given our model parameters, then maximize parameters given our expectations
                Expect(siteCounts, srcParameters, isError, config.BinomialTestThreshold, config.MinErrorRate);
                Maximize(siteCounts, isError, dstParameters, tmp, config.MinCoverage);

                // check if the parameters are sufficiently close to one another
                if (IsConverged(parameters, config.ConvergenceCriterion))
                {
                    converged = true;
                    finalParameters = dstParameters;
                    break;
                }
            }

            if (!converged)
                throw new ApplicationException("Variant detection algorithm could not converge");

            // predict using final parameter set (likely overkill)
            Expect(siteCounts, finalParameters, isError, config.BinomialTestThreshold, config.MinErrorRate);

        
            for (var from = 0; from < nBases; from++)
            {
                for (var to = 0; to < nBases; to++)
                {
                    if (from == to)
                        continue;
                     }
            }

            return ComputeVariants(siteCounts, finalParameters, isError, template,
                                   templateStart, config);
        }
    }
}
