using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Diagnostics;

namespace Bio.Util.Logging
{
    public class OutputInformation
    {

        private const double KB = 1024;
        private const double MB = KB * KB;
        private const double GB = MB * KB;

        /// <summary>
        /// Formats the specified memory in bytes to appropriate string.
        /// for example, 
        ///  if the value is less than one KB then it returns a string representing memory in bytes.
        ///  if the value is less than one MB then it returns a string representing memory in KB.
        ///  if the value is less than one GB then it returns a string representing memory in MB.
        ///  else it returns memory in GB.
        /// </summary>
        /// <param name="value">value in bytes</param>
        public static string FormatMemorySize(long value)
        {
            string result = null;
            if (value > GB)
            {
                result = (Math.Round(value / GB, 2)).ToString() + " GB";
            }
            else if (value > MB)
            {
                result = (Math.Round(value / MB, 2)).ToString() + " MB";
            }
            else if (value > KB)
            {
                result = (Math.Round(value / KB, 2)).ToString() + " KB";
            }
            else
            {
                result = value.ToString() + " Bytes";
            }

            return result;
        }
        public class MemorySize
        {
            public string PeakWorkingSet, TotalProcessorTime, PeakVirtualMemorySize64, PeakPagedMemorySize64, WorkingSet;
            public long NumericWorkingSet;
        }
        public static MemorySize GetMemoryUsage()
        {
              Process p = Process.GetCurrentProcess();
              return new MemorySize()
              {
                  PeakPagedMemorySize64 = FormatMemorySize(p.PeakPagedMemorySize64),
                  TotalProcessorTime = p.TotalProcessorTime.TotalSeconds.ToString(),
                  PeakVirtualMemorySize64 = FormatMemorySize(p.PeakVirtualMemorySize64),
                  PeakWorkingSet = FormatMemorySize(p.PeakWorkingSet64),
                  WorkingSet=FormatMemorySize(p.WorkingSet64),
                  NumericWorkingSet=p.WorkingSet64
              };
        }
            }
}
