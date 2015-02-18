using System;
using System.Collections.Generic;
using System.Linq;
using System.IO;
namespace ConstantModelOptimizer
{
    public static class EdnaToLisaConverter
    {
        /// <summary>
        /// Converter dictionary, based on the converter found in the debugger
        /// this is not guaranteed to be constant.
        /// </summary>
        public static Dictionary<int, char> upConverter = new Dictionary<int, char>() { {1,'T'}, {2,'G'}, {3,'A'}, {4,'C'}};
        public static StreamWriter sw = new StreamWriter("TemplatesAndReads.txt");
        public static void ConvertEdnaToLisa (string name, int[] template, int[] read)
        {
            var tpl = new string(template.Select(z=>upConverter[z]).ToArray());
            var rd = new string(read.Select(z=>upConverter[z]).ToArray());
            lock(sw) {
                sw.WriteLine(name);
                sw.WriteLine(tpl);
                sw.WriteLine(rd);
                sw.Flush();
            }

        }
    }
}

