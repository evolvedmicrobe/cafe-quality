﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Bio;
using Bio.Algorithms.Alignment.MultipleSequenceAlignment;
using Bio.SimilarityMatrices;

namespace VariantCaller
{
    /// <summary>
    /// This class represents a CCS MSA.  It is generated by adding several subreads to the 
    /// main read, and then examining 
    /// </summary>
    public class CCSMSA
    {
        /// <summary>
        /// Test currently fails, need to figure out what is wrong with aligner.
        /// </summary>
        public static void TestMuscleMultipleSequenceAlignment()
        {

            SimilarityMatrix similarityMatrix = new SimilarityMatrix(SimilarityMatrix.StandardSimilarityMatrix.AmbiguousDna);
            int gapOpenPenalty = -4;
            int gapExtendPenalty = -1;
            int kmerLength = 3;

            ISequence seqA = new Sequence(Alphabets.DNA, "GGGAAAAATCAGATT");
            ISequence seqB = new Sequence(Alphabets.DNA, "GGGAATCAAAATCAG");
            ISequence seqC = new Sequence(Alphabets.DNA, "GGGACAAAATCAG");
            List<ISequence> sequences = new List<ISequence>();
            sequences.Add(seqA);
            sequences.Add(seqB);
            sequences.Add(seqC);

            DistanceFunctionTypes distanceFunctionName = DistanceFunctionTypes.EuclideanDistance;
            UpdateDistanceMethodsTypes hierarchicalClusteringMethodName = UpdateDistanceMethodsTypes.Average;
            ProfileAlignerNames profileAlignerName = ProfileAlignerNames.NeedlemanWunschProfileAligner;
            ProfileScoreFunctionNames profileProfileFunctionName = ProfileScoreFunctionNames.WeightedInnerProduct;

            PAMSAMMultipleSequenceAligner msa = new PAMSAMMultipleSequenceAligner
                (sequences, kmerLength, distanceFunctionName, hierarchicalClusteringMethodName,
                profileAlignerName, profileProfileFunctionName, similarityMatrix, gapOpenPenalty, gapExtendPenalty,1,1);
//                Environment.ProcessorCount * 2, Environment.ProcessorCount);

            Console.WriteLine("Aligned sequences in stage 1: {0}", msa.AlignmentScoreA);
            for (int i = 0; i < msa.AlignedSequencesA.Count; ++i)
            {
                Console.WriteLine(new string(msa.AlignedSequencesA[i].Select(a => (char)a).ToArray()));
            }

            Console.WriteLine("Aligned sequences in stage 3: {0}", msa.AlignmentScoreC);
            for (int i = 0; i < msa.AlignedSequencesC.Count; ++i)
            {
                Console.WriteLine(new string(msa.AlignedSequencesC[i].Select(a => (char)a).ToArray()));
            }
            Console.WriteLine("Aligned sequences final: {0}", msa.AlignmentScore);

            for (int i = 0; i < msa.AlignedSequences.Count; ++i)
            {
                Console.WriteLine(new string(msa.AlignedSequences[i].Select(a => (char)a).ToArray()));
            }

            // Test case 2
            Console.WriteLine("Example 2");
            sequences = new List<ISequence>();
            sequences.Add(new Sequence(Alphabets.DNA, "GGGAAAAATCAGATT"));
            sequences.Add(new Sequence(Alphabets.DNA, "GGGAATCAAAATCAG"));
            sequences.Add(new Sequence(Alphabets.DNA, "GGGAATCAAAATCAG"));
            sequences.Add(new Sequence(Alphabets.DNA, "GGGAAATCG"));
            sequences.Add(new Sequence(Alphabets.DNA, "GGGAATCAATCAG"));
            sequences.Add(new Sequence(Alphabets.DNA, "GGGAATCTTATCAG"));
            sequences.Add(new Sequence(Alphabets.DNA, "GGGACAAAATCAG"));
            sequences.Add(new Sequence(Alphabets.DNA, "GGGAAAAATCAGATT"));
            sequences.Add(new Sequence(Alphabets.DNA, "GGGAATCAAAATCAG"));
            sequences.Add(new Sequence(Alphabets.DNA, "GGGACAAAATCAG"));


            msa = new PAMSAMMultipleSequenceAligner
                (sequences, kmerLength, distanceFunctionName, hierarchicalClusteringMethodName,
                profileAlignerName, profileProfileFunctionName, similarityMatrix, gapOpenPenalty, gapExtendPenalty,
                Environment.ProcessorCount * 2, Environment.ProcessorCount);

            Console.WriteLine("Aligned sequences in stage 1: {0}", msa.AlignmentScoreA);
            for (int i = 0; i < msa.AlignedSequencesA.Count; ++i)
            {
                Console.WriteLine(new string(msa.AlignedSequencesA[i].Select(a => (char)a).ToArray()));
            }

            Console.WriteLine("Aligned sequences in stage 3: {0}", msa.AlignmentScoreC);
            for (int i = 0; i < msa.AlignedSequencesC.Count; ++i)
            {
                Console.WriteLine(new string(msa.AlignedSequencesC[i].Select(a => (char)a).ToArray()));
            }
            Console.WriteLine("Aligned sequences final: {0}", msa.AlignmentScore);
            for (int i = 0; i < msa.AlignedSequences.Count; ++i)
            {
                Console.WriteLine(new string(msa.AlignedSequences[i].Select(a => (char)a).ToArray()));
            }

            // Test case e
            Console.WriteLine("Example 2");
            sequences = new List<ISequence>();
            sequences.Add(new Sequence(Alphabets.DNA, "GGGAAAAATCAGATT"));
            sequences.Add(new Sequence(Alphabets.DNA, "GGGAATCAAAATCAG"));
            sequences.Add(new Sequence(Alphabets.DNA, "GGGAATCAAAATCAG"));
            sequences.Add(new Sequence(Alphabets.DNA, "GGGAAATCG"));
            sequences.Add(new Sequence(Alphabets.DNA, "GGGAATCAATCAG"));
            sequences.Add(new Sequence(Alphabets.DNA, "GGGACAAAATCAG"));
            sequences.Add(new Sequence(Alphabets.DNA, "GGGAATCTTATCAG"));


            msa = new PAMSAMMultipleSequenceAligner
                (sequences, kmerLength, distanceFunctionName, hierarchicalClusteringMethodName,
                profileAlignerName, profileProfileFunctionName, similarityMatrix, gapOpenPenalty, gapExtendPenalty,
                Environment.ProcessorCount * 2, Environment.ProcessorCount);

            Console.WriteLine("Aligned sequences in stage 1: {0}", msa.AlignmentScoreA);
            for (int i = 0; i < msa.AlignedSequencesA.Count; ++i)
            {
                Console.WriteLine(new string(msa.AlignedSequencesA[i].Select(a => (char)a).ToArray()));
            }

            Console.WriteLine("Aligned sequences in stage 3: {0}", msa.AlignmentScoreC);
            for (int i = 0; i < msa.AlignedSequencesC.Count; ++i)
            {
                Console.WriteLine(new string(msa.AlignedSequencesC[i].Select(a => (char)a).ToArray()));
            }
            Console.WriteLine("Aligned sequences final: {0}", msa.AlignmentScore);
            for (int i = 0; i < msa.AlignedSequences.Count; ++i)
            {
                Console.WriteLine(new string(msa.AlignedSequences[i].Select(a => (char)a).ToArray()));
            }

        }

    }
}
