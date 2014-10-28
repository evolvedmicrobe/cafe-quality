using System;
using System.Reflection;
using System.Linq;
using System.Collections.Generic;

namespace VariantCaller
{
    public struct FieldsOfInterest {
        public FieldInfo[] NonArrayFields;
        public FieldInfo[] ArrayFields;
        public FieldsOfInterest(FieldInfo[] nonarrayfields, FieldInfo[] arrayfields)
        {
            NonArrayFields = nonarrayfields;
            ArrayFields = arrayfields;
        }
    }

    public class OutputHelper
    {
        /// <summary>
        /// Covariates prefixes, each of these is calculated for an array type.
        /// </summary>
        public static readonly string[] Array_Covariate_Prefixes = new [] {"Max_","Min_","Mean_"};

        public static readonly Func<IList<System.ValueType>, double>[] Array_Funcs = new Func<IList<ValueType>, double>[] {
            x => x.Max(p => (double) p),
            y => y.Min(p => (double) p),
            z => z.Average(p => (double) p)
        };


        /// <summary>
        /// Output the FieldInfos that will be used if this object is queried.
        /// </summary>
        /// <returns>The helper.</returns>
        /// <param name="needsHeader">Needs header.</param>
        /// <param name="doMaxMinMean">If set to <c>true</c> do max minimum mean.</param>
        static FieldsOfInterest FindFields (Type needsHeader, bool doMaxMinMean = false)
        {
            FieldInfo[] outputValues;
            var outAttribute = typeof(OutputAttribute);
            var outArray = typeof(OutputArrayAttribute);

            var outputNonArray = needsHeader.GetFields()
                            .Where(prop => prop.IsDefined(outAttribute, true))
                            .ToArray();

            var outputArrays = needsHeader.GetFields()
                .Where(prop =>prop.IsDefined(outArray,true))
                .ToArray();
            return new FieldsOfInterest(outputNonArray,outputArrays);
        }
        /// <summary>
        /// Generates an array of header columns for a particular class with fields labelled as output attributes.
        /// </summary>
        /// <returns>The headers.</returns>
        /// <param name="needsHeader">Needs header.</param>
        public static List<string> GetHeaders(Type needsHeader)
        {
            var fields = FindFields (needsHeader);
            List<string> headers = new List<string> (fields.ArrayFields.Length + fields.NonArrayFields.Length * Array_Covariate_Prefixes.Length);
            foreach (var f in fields.NonArrayFields) {
                headers.Add (f.Name);
            }
            foreach (var f in fields.ArrayFields) {
                foreach (var s in Array_Covariate_Prefixes) {
                    headers.Add (s + f.Name);
                }
            }
            return headers;
        }
        /// <summary>
        /// Actually gets the data lines for a particular class or object.
        /// </summary>
        /// <returns>The data lines.</returns>
        /// <param name="data">Data.</param>
        /// <typeparam name="T">The 1st type parameter.</typeparam>
        public static List<string> CalculateDataLines<T>(T data)
        {
            var t = data.GetType ();
            var fields = FindFields (t);
            List<string> dataFields = new List<string> (fields.ArrayFields.Length + fields.NonArrayFields.Length * Array_Covariate_Prefixes.Length);

            foreach (var f in fields.NonArrayFields) {
                dataFields.Add (GetValueSafe(data,f));
            }

            foreach (var f in fields.ArrayFields) {
                foreach (var s in Array_Funcs) {
                    try 
                    {
                        var value = f.GetValue (data) as IList<ValueType>;
                        var res = s (value);
                        dataFields.Add (res.ToString());
                    }
                    catch {
                        dataFields.Add ("NaN");
                    }
                }
            }
            return dataFields;
        }

        static string GetValueSafe<T> (T toReportOn, FieldInfo FI)
        {
            try {
                return FI.GetValue (toReportOn).ToString ();
            } catch {
                return "NULL";
            }
        }
    }
}

