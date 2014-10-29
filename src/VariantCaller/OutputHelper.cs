using System;
using System.Reflection;
using System.Linq;
using System.Collections.Generic;
using System.Collections;

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

        public static readonly Func<byte[], double>[] Array_Funcs_Byte = new Func<byte[], double>[] {
            x => x.Max(),
            y => y.Min(),
            z => z.Average( x=> (double)x)
        };
        public static readonly Func<short[], double>[] Array_Funcs_Short = new Func<short[], double>[] {
            x => x.Max(),
            y => y.Min(),
            z => z.Average(x => (double)x)
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
                var value = f.GetValue (data);
                if (value != null) {
                    // TODO: Clean up this code repetition.
                    var as_byte = value as byte[];
                    var as_short = value as short[];
                    if (as_byte != null) {
                        foreach (var s in Array_Funcs_Byte) {
                            try {
                                var res = s (as_byte);
                                dataFields.Add (res.ToString ());
                            } catch {
                                dataFields.Add ("NaN");
                            }
                        }
                    } else if (as_short != null) {
                        foreach (var s in Array_Funcs_Short) {

                            try {
                                var res = s (as_short);
                                dataFields.Add (res.ToString ());
                            } catch {
                                dataFields.Add ("NaN");
                            }
                        }

                    } else {
                        throw new Bio.BioinformaticsException ("missing converter");
                    }
                } else {
                    for (int i = 0; i < Array_Funcs_Byte.Length; i++) {
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

