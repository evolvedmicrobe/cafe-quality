using System;

namespace ConstantModelOptimizer
{
    public class DynamicProgrammingMatrixPair
    {
        /// <summary>
        /// Read position row, columns template
        /// </summary>
        public double[][] forward;
        /// <summary>
        /// Read position row, columns template
        /// </summary>
        public double[][] reverse;

        public DynamicProgrammingMatrixPair (string read, string template)
        {
            forward = RectangularArrays.ReturnRectangularDoubleArray (read.Length + 1, template.Length + 1);
        }
 
    }
}

