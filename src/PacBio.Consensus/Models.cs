using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Xml.Linq;
using ConsensusCore;
using IniParser;
using PacBio.Utils;

namespace PacBio.Consensus
{

    /// <summary>
    /// Handle loading the Quiver parameters from the QuiverParameters.ini file
    /// </summary>
    public static class ParameterLoading
    {
        private static PacBioLogger logger = PacBioLogger.GetLogger("ParameterLoading");

        private static void Log(LogLevel level, string msg)
        {
            logger.Log(level, msg);
        }

        private static void Log(LogLevel level, string msg, params object[] args)
        {
            logger.Log(level, String.Format(msg, args));
        }


        public static string InternalParameterDirectory
        {
            get
            {
                var path = Path.GetDirectoryName(new Uri(Assembly.GetExecutingAssembly().CodeBase).LocalPath);
                return Path.Combine(path, "Chemistry");
            }
        }

        public static string DefaultConsensusParametersPath
        {
            get { return InternalParameterDirectory; }
        }

        public static float ModelFactor = 1.0f;

        public static ScorerConfig DefaultQuiver
        {
            get
            {
                //TODO: Implement
                return null;
            }
        }

        public static ScorerConfig DefaultCCS
        {
            get
            {
                //TODO: Implement
                return null;
            }
        }

        /// <summary>
        /// Model Types in order of preference (we are not chekcing the availability of the different QVs -- just assume everything is available)
        /// </summary>
        private static List<string> ModelTypePreferences = new List<string>
            {
                "AllQVsMergingByChannelModel", "AllQVsModel", "NoMergeQVModel",
                "NoQVsModel", "NoQVsMergingByChannelModel"
            };

        public class ParameterSet
        {
            public string SequencingChemistry;
            public string ModelType;
        }


        public class ParameterException : Exception
        {
            public ParameterException()
            {}

            public ParameterException(string message)
                : base(message)
            {}

            public ParameterException(string message, Exception inner)
                : base(message, inner)
            {}
        }
        #if FALSE
        /// <summary>
        /// Load the Quiver parameters available in parameterFile, and return a QuiverConfigTable of all included parameter sets
        /// </summary>
        public static ScorerConfig LoadParametersFromFile(string parameterFile, string chemistry = null, string model = null)
        {
            if (!File.Exists(parameterFile))
            {
                throw new ApplicationException(String.Format("Parameter file: '{0}' doesn't exist", parameterFile));
            }

            var paramSets = new Dictionary<string, Tuple<int,QvModelParams>>();
            var unspecified = String.IsNullOrEmpty(chemistry) && String.IsNullOrEmpty(model);

            foreach (var paramSet in ReadParameters(parameterFile))
            {
                var rank = ModelTypePreferences.IndexOf(paramSet.ModelType);

                if (rank < 0)
                {
                    Log(LogLevel.WARN, "Unknown model type: '{0}', skipping...", paramSet.ModelType);
                    continue;
                }

                // empty chemistry and model, so load everything into the QuiverConfigTable
                if (unspecified)
                {
                    var key = paramSet.SequencingChemistry.ToLower().Contains("unknown")
                              ? "*" : paramSet.SequencingChemistry;

                    if (!paramSets.ContainsKey(key) || paramSets[key].Item1 > rank)
                        paramSets[key] = Tuple.Create(rank, paramSet.ModelParams);
                }
                // otherwise if we've just the chemistry, then load just it
                else if (String.IsNullOrEmpty(model) && paramSet.SequencingChemistry == chemistry)
                {
                    // if we don't have a key yet or our rank is superior, replace
                    if (!paramSets.ContainsKey("*") || paramSets["*"].Item1 > rank)
                        paramSets["*"] = Tuple.Create(rank, paramSet.ModelParams);
                }
                // otherwise if we're super-specific, just load that particular model
                else if (paramSet.SequencingChemistry == chemistry && paramSet.ModelType == model)
                {
                    paramSets["*"] = Tuple.Create(0, paramSet.ModelParams);
                }
            }

            if (paramSets.Count < 1)
            {
                var chemModel = chemistry ?? "autodetection";
                if (!String.IsNullOrEmpty(model))
                    chemModel += "." + model;
                throw new ParameterException(String.Format("No valid sequencing chemisty parameters for '{0}'", chemModel));
            }

            var table = new QuiverConfigTable();

            foreach (var entry in paramSets)
            {
                var scoreDiff = 18;
                var fastScoreThreshold = -12.5f;

                // if things are unknown, then loosen it up a bit
                if (entry.Key == "*")
                {
                    scoreDiff = 24;
                    fastScoreThreshold = -50.0f;
                }

                var config = new QuiverConfig(entry.Value.Item2, (int) Move.ALL_MOVES, new BandingOptions(4, scoreDiff), fastScoreThreshold);
                table.Insert(entry.Key, config);
            }

            Log(LogLevel.INFO, "Used parameter file: '{0}'", parameterFile);

            return new ScorerConfig {
                Parameters = table,
                Algorithm = RecursionAlgo.Viterbi,
                HasChemistryOverride = !unspecified
            };
        }


        public static string SelectParameterFile(string directoryOrFile, string parameterFileName)
        {
            string fn;
            
            // Case 1: No directory provided -- use SEYMOUR_HOME, or internal parameters
            if (String.IsNullOrEmpty(directoryOrFile))
            {
                fn = Path.Combine(DefaultConsensusParametersPath, parameterFileName);
                return fn;
            }

            // Case 2: File is fully specified on cmd-line
            if (File.Exists(directoryOrFile))
            {
                return directoryOrFile;
            }

            // Case 3: the path doesn't include GenomicConsensus
            if (Directory.Exists(directoryOrFile) &&
                Directory.GetDirectories(directoryOrFile).Select(Path.GetFileName).Contains("GenomicConsensus"))
            {
                return Path.Combine(directoryOrFile, "GenomicConsensus", parameterFileName);
            }

            // Case 4: Path is specified, append requested filename.
            return Path.Combine(directoryOrFile, parameterFileName);
        }


        private static object readerLock = new object();

        private static List<ParameterSet> ReadParameters(string filename)
        {
            lock (readerLock)
            {
                using (var sr = new StreamReader(filename))
                {
                    var parser = new StreamIniDataParser();
                    parser.CommentDelimiter = '#';

                    var parsedData = parser.ReadData(sr);
                    return parsedData.Sections.Select(ParseParams).ToList();
                }
            }
        }

        // I don't understand why we need a readerlock,
        // but I'll assume the need to use it during write too (IniFileParser isn't threadsafe?)

        public static void SaveModelToFile(string sequencingChemistry, string modelType, QvModelParams model, string filename)
        {
            Func<float, string> format = x => String.Format("{0:G9}", x / ModelFactor);

            using (var sr = new System.IO.StreamWriter(filename))
            {
                var iniData = new IniData();
                var sectionName = sequencingChemistry + "." + modelType;
                iniData.Sections.AddSection(sectionName);
                var sectionData = iniData.Sections.GetSectionData(sectionName);

                sectionData.Comments.AddRange(new string[] {"", " " + sequencingChemistry + " chemistry", ""});

                if (ModelTypePreferences.Where(m => m == modelType).Any())
                {
                    sectionData.Keys.AddKey("Match", format(model.Match));
                    sectionData.Keys.AddKey("Mismatch", format(model.Mismatch));
                    sectionData.Keys.AddKey("MismatchS", format(model.MismatchS));
                    sectionData.Keys.AddKey("Branch", format(model.Branch));
                    sectionData.Keys.AddKey("BranchS", format(model.BranchS));
                    sectionData.Keys.AddKey("DeletionN", format(model.DeletionN));
                    sectionData.Keys.AddKey("DeletionWithTag", format(model.DeletionWithTag));
                    sectionData.Keys.AddKey("DeletionWithTagS", format(model.DeletionWithTagS));
                    sectionData.Keys.AddKey("Nce", format(model.Nce));
                    sectionData.Keys.AddKey("NceS", format(model.NceS));
                }
                else
                {
                    throw new ApplicationException(String.Format("Unrecognized model type: {0}", modelType));
                }

                if (modelType == "AllQVsModel" || modelType == "NoMergeQVModel" || modelType == "NoQVsModel")
                {
                    sectionData.Keys.AddKey("Merge", format(model.Merge_A()));
                    sectionData.Keys.AddKey("MergeS", format(model.MergeS_A()));
                }
                else if (modelType == "AllQVsMergingByChannelModel" || modelType == "NoQVsMergingByChannelModel")
                {
                    sectionData.Keys.AddKey("Merge_A", format(model.Merge_A()));
                    sectionData.Keys.AddKey("Merge_C", format(model.Merge_C()));
                    sectionData.Keys.AddKey("Merge_G", format(model.Merge_G()));
                    sectionData.Keys.AddKey("Merge_T", format(model.Merge_T()));
                    sectionData.Keys.AddKey("MergeS_A", format(model.MergeS_A()));
                    sectionData.Keys.AddKey("MergeS_C", format(model.MergeS_C()));
                    sectionData.Keys.AddKey("MergeS_G", format(model.MergeS_G()));
                    sectionData.Keys.AddKey("MergeS_T", format(model.MergeS_T()));
                }

//                if (model is KineticModelParams)
//                {
//                    var kineticModel = (KineticModelParams) model;
//                    sectionData.Keys.AddKey("Mixing", format(kineticModel.Mixing));
//                    sectionData.Keys.AddKey("Rate", format(kineticModel.Rate));
//                }

                var parser = new StreamIniDataParser();
                parser.CommentDelimiter = '#';
                parser.WriteData(sr, iniData);
            }
        }
       
        private static ParameterSet ParseParams(SectionData section)
        {
            var f = ModelFactor;

            var name = section.SectionName.Split('.');

            var sequencingChemistry = name[0];
            var modelType = name[1];

            Func<string, float> GetParam = paramName => Single.Parse(section.Keys[paramName]);

            QvModelParams p;

            if (modelType == "AllQVsModel" || modelType == "NoMergeQVModel" || modelType == "NoQVsModel")
            {
                p = new QvModelParams(
                    GetParam("Match")*f,
                    GetParam("Mismatch")*f,
                    GetParam("MismatchS")*f,
                    GetParam("Branch")*f,
                    GetParam("BranchS")*f,
                    GetParam("DeletionN")*f,
                    GetParam("DeletionWithTag")*f,
                    GetParam("DeletionWithTagS")*f,
                    GetParam("Nce")*f,
                    GetParam("NceS")*f,
                    GetParam("Merge")*f,
                    GetParam("MergeS")*f);
            }
            else if (modelType == "AllQVsMergingByChannelModel" || modelType == "NoQVsMergingByChannelModel")
            {
                p = new QvModelParams(
                    GetParam("Match")*f,
                    GetParam("Mismatch")*f,
                    GetParam("MismatchS")*f,
                    GetParam("Branch")*f,
                    GetParam("BranchS")*f,
                    GetParam("DeletionN")*f,
                    GetParam("DeletionWithTag")*f,
                    GetParam("DeletionWithTagS")*f,
                    GetParam("Nce")*f,
                    GetParam("NceS")*f,
                    GetParam("Merge_A")*f,
                    GetParam("Merge_C")*f,
                    GetParam("Merge_G")*f,
                    GetParam("Merge_T")*f,
                    GetParam("MergeS_A")*f,
                    GetParam("MergeS_C")*f,
                    GetParam("MergeS_G")*f,
                    GetParam("MergeS_T")*f);
            }
            else
            {
                throw new ApplicationException(String.Format("Unrecognized model type: {0}", modelType));
            }

            return new ParameterSet
                {
                    SequencingChemistry = sequencingChemistry,
                    ModelType = modelType,
                    ModelParams = p
                };
        }
        #endif
    }

    public class ChemistryMapping
    {
        public static string GetMappingForMovie(string mappingFile, string movieName)
        {
            var mappings = GetChemistryMapping(mappingFile);
            return mappings[movieName];
        }

        /// <summary>
        /// Load the chemistry mappings from the chemistry xml file
        /// </summary>
        public static Dictionary<string, string> GetChemistryMapping(string mappingFile)
        {
            var mappingDict = new Dictionary<string, string>();
            var doc = XDocument.Load(mappingFile);
            doc.Descendants("Mapping");

            foreach (var n in doc.Descendants("Mapping"))
            {
                var movie = n.Descendants("Movie").FirstOrDefault().Value;
                var sequencingChemistry = n.Descendants("SequencingChemistry").FirstOrDefault().Value;
                mappingDict[movie] = sequencingChemistry;
            }

            return mappingDict;
        }


        public static string GetSequencingChemistryForMovie(string movieFile, string chemistry)
        {
            var sequencingChemistry = "unknown";

            // figure out what the current movie name is.
            // movieName may have up to 3 extensions
            var movieName = movieFile;
            movieName = Path.GetFileNameWithoutExtension(movieName);
            movieName = Path.GetFileNameWithoutExtension(movieName);
            movieName = Path.GetFileNameWithoutExtension(movieName);

            if (File.Exists(chemistry))
            {
                // if we have a chemistry_mapping.xml, load it.
                var mapping = ChemistryMapping.GetChemistryMapping(chemistry);

                if (!mapping.ContainsKey(movieName))
                {
                    throw new ApplicationException(String.Format("Chemistry mapping File: '{0}' doesn't contain mapping for movie: '{1}'", chemistry, movieFile));
                }

                sequencingChemistry = mapping[movieName];
            }
            else
            {
                // if the Chemistry argument isn't a file, interpret it as a SequencingChemistry
                sequencingChemistry = chemistry;
            }

            return sequencingChemistry;
        }
    }
}
