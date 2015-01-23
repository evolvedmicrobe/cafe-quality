﻿using System;
using System.Collections.Generic;
using System.Linq;
using ConsensusCore;
using PacBio.IO;
using PacBio.Utils;
using PacBio.Align;

namespace PacBio.Consensus
{
    /// <summary>
    /// Temporary bandaid class to determine if we can improve homopolymers by proposing mutations in homopolymers >4 and 
    /// scoring them with the C2 chemistry. 
    /// </summary>
    public static class HomopolymerPolisher
    {
        public static int TotalMutationsAccepted = 0;
        //Early work showed these changes work the best on data with an SNR below 9
        public const double MeanSNRCutoff = 9.0;
        public static Tuple<TrialTemplate, List<MutationScore>> PolishHomopolymers(TrialTemplate tpl,
                                                            MultiReadMutationScorer oldScorer,
                                                            IZmwBases bases,
                                                            List<MutationScore> allScores)
        {
            var meanSNR = bases.Metrics.HQRegionSNR.Average ();
            // Don't try to fix with SNR > 9, you introduce more errors.
            if (meanSNR > MeanSNRCutoff) {
                return new Tuple<TrialTemplate, List<MutationScore>>(tpl, allScores);
            }
            var scConfig = ParameterLoading.C2Parameters;
            scConfig.Algorithm = RecursionAlgo.Prob;
            var scorer = new MultiReadMutationScorer(oldScorer.OriginalRegions, bases, tpl, scConfig);

            Func<Mutation, MutationScore> scoreMutation =
                m => new MutationScore { Score = scorer.ScoreMutation(m), Mutation = m, Exists = true };
            double score;
            int mutationSpacing = 1;
            // This corresponds to a probability of ~51% - to avoid flip-flopping between opposite mutations
            float minScore = 0.35f;
            Func<IEnumerable<Mutation>, List<Mutation>> screenMutations = mutationsToTry => FindConsensus.FindMutations(mutationsToTry, scoreMutation, out score, mutationSpacing, minScore);
            List<Mutation> mutsToTry =  GenerateLongHomopolymerMutations(tpl).ToList();
            mutsToTry.Reverse();
            var accepted = screenMutations(mutsToTry);
            List<MutationScore> newScores = new List<MutationScore>();
            foreach (var m in accepted)
            {
                for ( int i=0; i<allScores.Count; i++)
                {
                    var mut = allScores[i];
                    if (mut.Mutation.TemplatePosition > m.TemplatePosition)
                    {
                        mut.Mutation.TemplatePosition++;
                        allScores[i] = mut;
                    }
                }
            }
            for (int i=0; i< accepted.Count; i++)
            {
                System.Threading.Interlocked.Increment (ref TotalMutationsAccepted);
                if (i % 100 == 0) {
                    Console.WriteLine ("Changed: " + TotalMutationsAccepted + " mutations");
                }
            }
            scorer.ApplyMutations(accepted);
            return new Tuple<TrialTemplate, List<MutationScore>>(scorer.Template, allScores);
        }
       static IEnumerable<Mutation> GenerateLongHomopolymerMutations(TrialTemplate tpl, int minLength = 2)
        {
            // Attempt to insert or mismatch every base at every positions.
            // Don't mutate anything in the know template adapter region
            var seq = tpl.GetSequence(Strand.Forward);

            // Start is the first base past the adapter
            var start = Math.Max(1, tpl.StartAdapterBases);

            // End is the fist base of the end adapter
            var end = seq.Length - Math.Max(1, tpl.EndAdapterBases);

            for (int i = start; i <= end; i++)
            {

                // homopolyerStart is true if we are on the first base of a hp
                // We will only do a matching insertion or a deletion if we're a 
                // the start of a homopolyer, to prevent retesting the same template
                bool homopolyerStart = (i > 1 && seq[i - 1] == seq[i]);
                if (homopolyerStart)
                {
                    int len = 1;
                    int start_hp = i;
                    char bp = seq[i];
                    while (i <= end && seq[i] == bp)
                    {
                        len++;
                        i++;
                    }
                    if (len >= minLength && (bp == 'G' || bp == 'C'))
                    {
                        yield return new Mutation() { Base = bp, TemplatePosition = start_hp, Type = MutationType.INSERTION };
                    }
                    i--;
                }
            }     
        }


       
    }
}