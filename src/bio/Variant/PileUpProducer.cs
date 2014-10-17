using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Bio;
using Bio.IO.SAM;
using System.Diagnostics;

namespace Bio.Variant
{
    /// <summary>
    /// A class the produces read pile ups (columns of multiple sequence alignments) from a stream of 
    /// SAM aligned quality sequences.
    /// </summary>
	public static class PileUpProducer
    {        
        /// <summary>
        /// Takes a sorted list of SAMAligned sequences and converts them in to read pile-ups. One pile-up for genomic position or alignment.
        /// </summary>
        /// <param name="sequences"></param>
        /// <returns></returns>
		public static IEnumerable<PileUp> CreatePileupFromReads(IEnumerable<CompactSAMSequence> sequences)
		{
			LinkedList<PileUp> pileupsToEmit = null;
			string cur_ref, last_ref = null;

			// Variable for catching errors, can be removed later after once no problems show up.
			int lastEmittedPosition = -1;
			foreach (var seq in sequences) {
				if (!validateSequence (seq)) {
					continue;
				}
				;
				cur_ref = seq.RName;                

				// initalize if first sequence or a new reference
				if (pileupsToEmit == null || last_ref != cur_ref) {
					if (pileupsToEmit != null) {
						foreach (var pu in pileupsToEmit) {
							Debug.Assert (pu.Position >= lastEmittedPosition, "Emitting pile-up out of order");
							lastEmittedPosition = pu.Position;
							yield return pu;
						}
					}
					last_ref = cur_ref;
					lastEmittedPosition = -1;
					pileupsToEmit = new LinkedList<PileUp> ();
				}

				// Get read data to add to the pile-up.
				var bases = getBasesForSequence (seq);

				// Advance the pile-up until we start to overlap.
				while (pileupsToEmit.Count > 0 && pileupsToEmit.First.Value.Position < bases [0].Position) {
					var pu = pileupsToEmit.First.Value;
					Debug.Assert (pu.Position >= lastEmittedPosition, "Emitting pile-up out of order");
					lastEmittedPosition = pu.Position;
					yield return pu;
					pileupsToEmit.RemoveFirst ();
				}
                
				// Now add each of this read's bases to the pile-up
				var cur_pileup = pileupsToEmit.First;
				foreach (var bp in bases) {   
					if (cur_pileup == null) {
						// If null, then this base and all following bases are new additions.
						pileupsToEmit.AddLast (new PileUp (cur_ref, bp));                         
					} else {
						// Merge the two bases and the pile-up
						while (cur_pileup != null) {
							if (cur_pileup.Value.Position == bp.Position && cur_pileup.Value.InsertionOffSet == bp.InsertionOffSet) {
								cur_pileup.Value.Bases.Add (bp.BaseWithQuality);
								cur_pileup = cur_pileup.Next;
								break;
							} else if (bp.Position < cur_pileup.Value.Position ||
							                              (bp.Position == cur_pileup.Value.Position && bp.InsertionOffSet < cur_pileup.Value.InsertionOffSet)) {
								pileupsToEmit.AddBefore (cur_pileup, new PileUp (cur_ref, bp));
								break;
							} else {
								cur_pileup = cur_pileup.Next;
								if (cur_pileup == null) {
									pileupsToEmit.AddLast (new PileUp (cur_ref, bp));
								}
								break;
							}
						}
					}
				}
			}

			// Finally emit anything we haven't yet.
			if (pileupsToEmit != null) {
				foreach (var pu in pileupsToEmit) {
					Debug.Assert (pu.Position >= lastEmittedPosition, "Emitting pile-up out of order");
					lastEmittedPosition = pu.Position;
					yield return pu;
				}
			}
		}

        /// <summary>
        /// Method throws an exception if sequence violates any assumption made by this class anywhere.
        /// Avoids, separate checks within each method.
        /// </summary>
        /// <param name="seq"></param>
		private static bool validateSequence(CompactSAMSequence seq)
        {
            if (seq == null) {
                throw new ArgumentNullException("seq");
            }
            if (String.IsNullOrEmpty(seq.RName) || 
				seq.RefEndPos <= seq.Pos || 
                String.IsNullOrEmpty(seq.CIGAR) || 
				seq.CIGAR == "*" )
            {
				return false;
				//throw new ArgumentException("Tried to build a pileup with an invalid sequence.  Sequence was:\n"+
				//    seq.ToString());
            }
			return true;
        }

		        /// <summary>
        /// Turn a SAMAlignedSequence into a list of BaseAndQualityAndPosition objects,
        /// useful for adding to a pile-up.
        /// </summary>
        /// <param name="seq"></param>
        /// <returns></returns>
		static List<BaseAndQualityAndPosition> getBasesForSequence(CompactSAMSequence seq)
        {
            List<BaseAndQualityAndPosition> toReturn = new List<BaseAndQualityAndPosition>(seq.RefEndPos - seq.Pos + 10);
            // Decode the cigar string into operations.
            // TODO: This code is duplicated in many places
            string CIGAR = seq.CIGAR;
            List<KeyValuePair<char, int>> charsAndPositions = new List<KeyValuePair<char, int>>();
            for (int i = 0; i < CIGAR.Length; i++)
            {
                char ch = CIGAR[i];
                if (Char.IsDigit(ch))
                {
                    continue;
                }
                charsAndPositions.Add(new KeyValuePair<char, int>(ch, i));
            }

            // Get sequence bases and error probabilities
			var seq_phred_scores = seq.GetPhredQualityScores();
            var seq_bases = seq.ToArray();
            // Use the cigar operations to emit bases.
            int curRef = seq.Pos;
            int curQuery = 0;
            for (int i = 0; i < charsAndPositions.Count; i++)
            {
                // Parse the current cigar operation
                char ch = charsAndPositions[i].Key;
                int cig_start = i==0 ? 0 : charsAndPositions[i - 1].Value + 1;
                int cig_end = charsAndPositions[i].Value - cig_start;
                int cig_len = int.Parse(CIGAR.Substring(cig_start, cig_end));
                // Emit or advance based on cigar operation.
                switch (ch)
                {
                    case 'P': //padding (Silent deltions from padded reference)
                    case 'N': //skipped region from reference
                        throw new Exception("Pile up methods not built to handle reference clipping (Cigar P or N) yet.");
                    case 'M': //match or mismatch
                    case '=': //match
                    case 'X': //mismatch
                        for (int k = 0; k < cig_len; k++)
                        {                            
							var bqp= new BaseAndQualityAndPosition(curRef,0, new BaseAndQuality(seq_bases[curQuery], (byte) seq_phred_scores[curQuery]));
                            toReturn.Add(bqp);
                            curQuery++;
                            curRef++;
                        }
                        break;
                    case 'I'://insertion to the reference
                        for (int k = 0; k < cig_len; k++)
                        {                            
						var bqp =  new BaseAndQualityAndPosition(curRef,k, new BaseAndQuality(seq_bases[curQuery], (byte) seq_phred_scores[curQuery]));
                            toReturn.Add(bqp);
                            curQuery++;
                        }
                        break;
                    case 'D'://Deletion from the reference
                        for (int k = 0; k < cig_len; k++)
                        {                            
						var bqp = new BaseAndQualityAndPosition(curRef,k, new BaseAndQuality((byte)'-', byte.MinValue));
                            toReturn.Add(bqp);
                            curRef++;
                        }
                        break;
                    case 'S': //soft clipped
                        curQuery += cig_len;
                        break;
                    case 'H'://had clipped
                        break;
                    default:
                        throw new FormatException("Unexpected SAM Cigar element found " + ch.ToString());
                }                
            }
            return toReturn;
        }
        
    }
}
