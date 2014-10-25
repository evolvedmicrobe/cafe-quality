using System;

namespace VariantCaller
{
    /// <summary>
    /// Used for decorating scalar quantities that we are outputting to a file.
    /// e.g. RecordCount, TotalIterations, stuff like that.
    /// </summary>
    public class OutputAttribute : System.Attribute
    {
       
    }
    /// <summary>
    /// Used to decorate fields that are amenable to calculating max, min, mean, etc.
    /// values on.
    /// </summary>
    public class OutputArrayAttribute : System.Attribute {}
}

