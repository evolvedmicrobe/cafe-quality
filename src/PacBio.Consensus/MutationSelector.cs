using System;
using System.Collections.Generic;
using System.Linq;
using PacBio.Utils;


namespace PacBio.Consensus
{
    /// <summary>
    /// We have a list of mutations. Each mutation, is made up of a score, and a template position.
    /// This class selects the set of mutations that are all least n bases separated that have the highest
    /// total score.
    /// </summary>
    public class SpacedSelector
    {
        /// <summary>
        /// A convenience function for finding the best set of mutations that are spaced by at least <code>spacing</code>
        /// template positions
        /// </summary>
        public static List<MutationScore> BestMutations(List<MutationScore> mutations, int spacing)
        {
            return BestItems(mutations, ms => ms.Mutation.TemplatePosition, ms => ms.Score, spacing);
        }

        /// <summary>
        /// Internal class to aid in the BestItems method below.
        /// </summary>
        internal struct ItemToSort<T>
        {
            internal readonly T Item;
            internal readonly int Position;
            internal readonly double Score;
            public ItemToSort(T item, int position, double score)
            {
                Item = item;
                Score = score;
                Position = position;
            }
        }

        /// <summary>
        /// We are given a list of items. For each item we can compute a position and a score.  This method finds the set
        /// of items with the maximum total score such that are all at least <code>spacing</code> units away from one another
        /// 
        /// Important!: Positions and scores for items are assumed to be unchanged inside the method.
        /// </summary>
        /// <typeparam name="T">The item type</typeparam>
        /// <param name="items">The set of items</param>
        /// <param name="positionFunc">A function that returns the position of an item</param>
        /// <param name="scoreFunc">A function that returns the score of an item</param>
        /// <param name="minSpacing">The minimum spacing required between items in the returned set</param>
        /// <returns>The best set of items</returns>
        public static List<T> BestItems<T>(IEnumerable<T> items, Func<T, int> positionFunc, Func<T, double> scoreFunc, int minSpacing)
        {
            var itemArray = items.Select(x => new ItemToSort<T>(x, positionFunc(x), scoreFunc(x))).ToList();
            itemArray.Sort((x, y) => x.Position.CompareTo(y.Position));
            
            if (itemArray.Count == 0)
                return new List<T>();
            
            // The best score achieved by including item i
            double[] score = new double[itemArray.Count];
            // The previous included item if item i is included
            int[] prevItem = new int[itemArray.Count];

            // The first item has only it's own score, and no predecessor
            score[0] = itemArray[0].Score;
            prevItem[0] = -1;

            // For each subsequent item, figure out what the best preceding item would be
            for(int i = 1; i < itemArray.Count; i++)
            {
                // The starting point is to have no predecessor, and only item[i]'s score
                score[i] = itemArray[i].Score;
                prevItem[i] = -1;

                // Go through previous items, and find the one that satisfies the spacing constraint,
                // and has the highest score.
                for (int j = i - 1; j >= 0; j--)
                {
                    if (itemArray[i].Position - itemArray[j].Position > minSpacing)
                    {
                        var newScore = score[j] + itemArray[i].Score;
                        if (newScore > score[i])
                        {
                            score[i] = newScore;
                            prevItem[i] = j;
                        }
                    }
                }
            }

            // Find the cell with the highest total score, then trace back through the included items.
            var end = score.IMax();
            var final = new List<T>();

            while (end >= 0)
            {
                final.Insert(0, itemArray[end].Item);
                end = prevItem[end];
            }

            return final;
        }
    }
}
