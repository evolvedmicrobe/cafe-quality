using System;
using Bio;
using PacBio.Consensus;
using PacBio.Data;
using ConsensusCore;

namespace VariantCaller
{

    enum ScoringResult { Success, NoOverlapOrAlignment, BadData, MysteriousFail }
    enum ScoringMethod { SumProduct, Viterbi }
    /// <summary>
    /// This class represents a section of a template.  It is designed 
    /// to be useful for examining how different reads align to different templates.
    /// </summary>
    public class TemplateRegion
    {
        /// <summary>
        /// The full template.
        /// </summary>
        public Reference RefSeq;
        /// <summary>
        /// The start of this region on the original template. 0 indexed.
        /// </summary>
        public int StartRegion;
        /// <summary>
        /// The end of the region on the template, 0 indexed, not inclusive (so +1 the last position).
        /// </summary>
        public int EndRegion;
        /// <summary>
        /// The sequence of just the region
        /// </summary>
        public Sequence Haplotype;


        private static QuiverConfig q_config;
        public static TemplateRegion()
        {
            // This obviously can't stand.
            var ccsFileLocation = @"/Users/nigel/git/cafe-quality/data/CCSParameters.ini";
            var scoringConfig = ParameterLoading.LoadParametersFromFile(ccsFileLocation,"P6-C4");
            q_config = scoringConfig.Parameters.At("P6-C4");
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="VariantCaller.TemplateRegion"/> class.
        /// </summary>
        /// <param name="refseq">Refseq.</param>
        /// <param name="start">Start. 0 - indexed</param>
        /// <param name="end">End. NON-INCLUSIVE, that is 1 based the last position.  0-indexed.</param>
        public TemplateRegion (Reference refseq, int start, int end)
        {
            RefSeq = refseq;
            Haplotype = refseq.RefSeq.GetSubSequence(start, end - start) as Sequence;
            StartRegion = start;
            EndRegion = end;
            ReadScorer rs = new ReadScorer (q_config);
              
        }

        public ScoringResult ScoreReadAgainstTemplateRegion(CCSRead read, ScoringMethod method, out float logScore)
        {
            if (read.AssignedReference != RefSeq) {
                return ScoringResult.NoOverlapOrAlignment;
            } else if (read.ZMW == null) {
                return ScoringResult.BadData;
            }

            // First we want to align the sequence and see if we have overlap with reference.
            var cdata = read.ZMW.FullRead;

            //Val

            read.ZMW.FullRead;

            var score = new SparseSseQvMultiReadScorer (q_config);
            //score.Score(



            //using (var qsf = ConsensusCoreWrap.MakeQvSequenceFeatures(r.Start, r.End - r.Start, bases))
            //using (var read = new Read(qsf, name, chem))


               
        }

        private void Validate


    }
}

