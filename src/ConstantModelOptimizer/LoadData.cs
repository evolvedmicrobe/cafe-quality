using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;

namespace ConstantModelOptimizer
{
    /// <summary>
    /// This class is designed to grab ALL4MERS data to fit it for a particular SNR ranges
    /// </summary>
    public class LoadData
    {
        //public static Dictionary<int, ZmwInfo> LoadSNRs ()
        public static IEnumerable<ZmwInfo> LoadSNRs ()

        {
            return File.ReadLines ("master_5ba5286_combined_reads.csv").Skip (1).
                Select (z => new ZmwInfo (z)).
                Where (z => z.Reference == "ALL4MER.V2.01");//
                //.ToDictionary(z=>z.HoleNumber, y=>y);
            //return data;
        }
        public static IEnumerable<ReadTemplateInfo> LoadSampleData()
        {
            return File.ReadLines ("TemplateReadPairs.csv").Skip(1).Select (z => new ReadTemplateInfo (z));
        }

    }

    public class ReadTemplateInfo {
        public string read;
        public string template;
        public int subread;
        public int Hole;
        public ReadTemplateInfo(string line)
        {
            var sp = line.Trim().Split (',');
            Hole = Convert.ToInt32 (sp [0]);
            subread = Convert.ToInt32 (sp [1]);
            template = sp [2];
            read = sp [3];
        }

    }

    public class ZmwInfo
    {
        public string Reference;
        public int HoleNumber;
        public double SnrT, SnrG, SnrA,SnrC;
        public ZmwInfo(string line)
        {
            var sp = line.Split (',');
            Reference = sp [2];
            HoleNumber = Convert.ToInt32 (sp [9]);
            var doubles = sp.Skip (11).Take (4).Select (p => Convert.ToDouble (p)).ToArray ();
            SnrT = doubles [0];
            SnrG = doubles [1];
            SnrA = doubles [2];
            SnrC = doubles [3];
        }
    }
}

