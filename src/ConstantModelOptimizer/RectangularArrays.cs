using System;
using System.Collections;
internal static partial class RectangularArrays
{
    internal static double[][] ReturnRectangularDoubleArray(int Size1, int Size2)
    {

        double[][] Array;
        Array = new double[Size1][];
        for (int Array1 = 0; Array1 < Size1; Array1++)
        {
            Array[Array1] = new double[Size2];
        }
            
        
        return Array;
    }
    internal static double[][] ReturnRectangularDoubleFilledBad(int Size1, int Size2)
    {

        double[][] Array;
        Array = new double[Size1][];
        for (int Array1 = 0; Array1 < Size1; Array1++)
        {
            var arr = new double[Size2];
            fill (arr, Double.NegativeInfinity);
            Array[Array1] = arr;
        }
        return Array;
    }

    internal static void fill(double[] arr, double val)
    {
        for(int i=0;i<arr.Length;i++)
        {
            arr[i]=val;
        }
    }
}