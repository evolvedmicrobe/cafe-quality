using System;

namespace ConstantModelOptimizer
{
    public class DynamicProgrammingMatrixPair
    {
        /// <summary>
        /// Read position row, columns template
        /// </summary>
        public LatentStates[][] Forward;
        /// <summary>
        /// Read position row, columns template
        /// </summary>
        public LatentStates[][] Reverse;

        public DynamicProgrammingMatrixPair (string read, string template)
        {
            Forward = RectangularArrays.ReturnRectangularLatentStateArray (read.Length, template.Length);
            Reverse = RectangularArrays.ReturnRectangularLatentStateArray (read.Length, template.Length);

        }
        /// <summary>
        /// Empty the arrays and fill with negative infinity
        /// </summary>
        public void Clear()
        {
            RectangularArrays.ClearLatentArray (Forward);
            RectangularArrays.ClearLatentArray (Reverse);

        }
 
    }
}

