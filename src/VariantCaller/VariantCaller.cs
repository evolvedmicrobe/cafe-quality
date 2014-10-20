using System;
using System.Linq;
using System.Collections.Generic;
using Bio;
using Bio.Algorithms.Alignment;

namespace VariantCaller
{
	public static class VariantCaller
	{
//		public IList<Variant> CallVariants(PairwiseAlignedSequence alignment)
//		{
//			List<Variant> variants = new List<Variant>();
//			var reference = alignment.FirstSequence;
//			var query = alignment.SecondSequence;
//			int index;
//
//			for (index = 0; index < reference.Count; index++)
//			{
//				byte referenceCharacter = reference[index];
//				byte queryCharacter = query[index];
//
//				if (DnaAlphabet.Instance.Gap != referenceCharacter
//					&& DnaAlphabet.Instance.Gap != queryCharacter)
//				{
//					if (referenceCharacter != queryCharacter) {
//						alignment.FirstOffset
//					}
//				}
//				else
//				{		
//					int gapCount;
//					if (DnaAlphabet.Instance.Gap == referenceCharacter)
//					{
//						gapCount = FindExtensionLength(referenceSequence, index);
//					}
//					else
//					{
//						gapCount = FindExtensionLength(querySequence, index);
//					}
//					score += GapOpenCost + (gapCount * GapExtensionCost);
//
//					// move the index pointer to end of extension
//					index = index + gapCount - 1;
//				}
//			}
//			return score;
//		}
	}
}

