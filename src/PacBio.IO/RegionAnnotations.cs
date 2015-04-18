#region Copyright (c) 2010, Pacific Biosciences of California, Inc.
//
// All rights reserved.
//
// THIS SOFTWARE CONSTITUTES AND EMBODIES PACIFIC BIOSCIENCES’ CONFIDENTIAL
// AND PROPRIETARY INFORMATION.
//
// Disclosure, redistribution and use of this software is subject to the
// terms and conditions of the applicable written agreement(s) between you
// and Pacific Biosciences, where “you” refers to you or your company or
// organization, as applicable.  Any other disclosure, redistribution or
// use is prohibited.
//
// THIS SOFTWARE IS PROVIDED BY PACIFIC BIOSCIENCES AND ITS CONTRIBUTORS "AS
// IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
#endregion

using System;
using System.Collections.Generic;
using System.Linq;
using PacBio.IO.ICD;
using PacBio.HDF;
using PacBio.Utils;

namespace PacBio.IO
{

    public interface IZmwRegions
    {
        IList<RegionAnnotator.Region> Regions { get; }
        ISequencingZmw Zmw { get; }
    }
    
    /// <summary>
    /// Tools for reading and writing trace region annotations to the Region table in the bas.h5 file.
    /// The region table lives at /PulseData/Regions.
    /// </summary>
    [HdfIcdEntry(Path = "/PulseData")]
    public class RegionAnnotator
    {
        #region Interface Control Documentation

        public class Icd : HdfIcd<RegionAnnotator>
        {
            // The Regions table is fully documented here
            private static IEnumerable<HdfIcdEntry> Mine()
            {
                return new[]
                    {
                        new HdfIcdEntry
                            {
                                Path = "Regions",
                                Detail = "The Regions table, with columns identifed by ColumnNames attribute"
                            },
                        new HdfIcdEntry {Path = "Regions/RegionTypes", Detail = "Region type lookup table"},
                        new HdfIcdEntry {Path = "Regions/RegionDescriptions", Detail = "Region type description"},
                        new HdfIcdEntry
                            {
                                Path = "Regions/RegionSources",
                                Detail = "Origin or source of the region annotation"
                            },
                        new HdfIcdEntry
                            {
                                Path = "Regions/ColumnNames",
                                Detail = "Identification of Regions table columns"
                            },
                    };
            }

            public Icd(bool flatten = false)
                : base(Mine(), flatten)
            {
            }
        }

        #endregion

        public static RegionType AdapterRegionType
            = new RegionType
                {
                    Description = "Adapter Hit",
                    Source = "AdapterFinding",
                    Type = "Adapter"
                };

        public static RegionType InsertRegionType =
            new RegionType
                {
                    Description = "Insert Region",
                    Source = "AdapterFinding",
                    Type = "Insert"
                };

        public static RegionType HqRegionType =
            new RegionType
                {
                    Type = "HQRegion",
                    Description =
                        "High Quality bases region. Score is 1000 * predicted accuracy, where predicted accuary is 0 to 1.0",
                    Source = "PulseToBase Region classifer"
                };

        /// <summary>
        /// Create a Region table reader for the already-opened HDF5 file f
        /// </summary>
        public static Reader GetReader(IChunkFile f)
        {
            var annotationGroup = (IGroup) f.File.GetChild("PulseData");

            if (annotationGroup == null)
                return null;

            return new Reader(annotationGroup);
        }

        /// <summary>
        /// Describes a class of regions that may exist in the Region table
        /// </summary>
        public class RegionType : IEquatable<RegionType>
        {
            /// <summary>
            /// Region type name (e.g. Adapter, Insert, etc)
            /// </summary>
            public string Type;

            /// <summary>
            /// Description of the meaning of the region type
            /// </summary>
            public string Description;

            /// <summary>
            /// Description of the algorithm or component that generated this region type
            /// </summary>
            public string Source;

            public override int GetHashCode()
            {
                return String.Concat(Type, Description, Source).GetHashCode();
            }

            // Required to do the right lookup in the dictionary
            public bool Equals(RegionType other)
            {
                return GetHashCode() == other.GetHashCode();
            }
        }

        public static readonly string RegionDatasetName = "Regions";
        public static string TypesAttribute = "RegionTypes";
        public static string DescriptionAttribute = "RegionDescriptions";
        public static string SourceAttribute = "RegionSources";

        public static int HoleNumberCol = 0;
        public static int RegionTypeIdCol = 1;
        public static int RegionStartCol = 2;
        public static int RegionEndCol = 3;
        public static int ScoreCol = 4;

        /// <summary>
        /// Simple all-at-once region table reader
        /// </summary>
        public class Reader
        {
            public Reader(IGroup annotationGroup)
            {
                // Get the region dataset
                var dataset = (IDataset) annotationGroup.GetChild(RegionDatasetName);

                // Read the region types present
                var regionTypeNames = (string[]) dataset.GetAttribute(TypesAttribute).Read();
                var regionDescriptions = (string[]) dataset.GetAttribute(DescriptionAttribute).Read();
                var regionSources = (string[]) dataset.GetAttribute(SourceAttribute).Read();

                // Write copy the region table into our region structs
                regionTypes = regionTypeNames.Length.Fill(i =>
                                                          new RegionType
                                                              {
                                                                  Description = regionDescriptions[i],
                                                                  Source = regionSources[i],
                                                                  Type = regionTypeNames[i]
                                                              });

                var rawAnnotations = (int[,]) dataset.Read();
                Regions = MakeRegions(rawAnnotations);
            }

            private RegionType[] regionTypes;

            /// <summary>
            /// All the region annotation in a base file
            /// </summary>
            public Region[] Regions { get; private set; }

            private Region[] MakeRegions(int[,] ra)
            {
                return ra.GetLength(0).Fill(i =>
                                            new Region
                                                {
                                                    End = ra[i, RegionEndCol],
                                                    Start = ra[i, RegionStartCol],
                                                    Score = ra[i, ScoreCol],
                                                    HoleNumber = ra[i, HoleNumberCol],
                                                    Type = regionTypes[ra[i, RegionTypeIdCol]]
                                                });
            }
        }

        /// <summary>
        /// 
        /// </summary>
        public class Region
        {
            /// <summary>
            /// The Zmw HoleNumber is region annotation applies to
            /// </summary>
            public int HoleNumber;

            /// <summary>
            /// The first basecall in the region
            /// </summary>
            public int Start;

            /// <summary>
            /// The index after the last basecall in the region
            /// </summary>
            public int End;

            /// <summary>
            /// The score of the region
            /// </summary>
            public int Score;

            /// <summary>
            /// A reference to the RegionType descriptor of the region
            /// </summary>
            public RegionType Type;

            public int Length
            {
                get { return End - Start; }
            }

        }


        /// <summary>
        /// Write Region annotations to a Region dataset in an HDF5 file.
        /// </summary>
        public class Writer
        {
            private IDataset regionDataset;
            private Dictionary<RegionType, int> regionTypeCodes = new Dictionary<RegionType, int>();
            private List<RegionType> regionTypeList = new List<RegionType>();

            public Writer(IGroup targetGroup)
            {
                // Set up the Region dataset, and the attributes describing the region types in use
                regionDataset = (IDataset) targetGroup.GetChild(RegionDatasetName);


                if (regionDataset == null)
                {
                    // If it doesn't already exist, setup the new dataset
                    var dataType = targetGroup.File.CreateDatatype(typeof (int));
                    var dataSpace = targetGroup.File.CreateDataspace(new long[] {0, 5}, new long[] {-1, 5});
                    regionDataset = targetGroup.CreateDataset(RegionDatasetName, dataType, dataSpace);

                    // Attribute describing the meaning of the Region table columns
                    regionDataset.InsertAttribute("ColumnNames",
                                                  new[]
                                                      {
                                                          "HoleNumber", "Region type index", "Region start in bases",
                                                          "Region end in bases", "Region score"
                                                      });
                }
                else
                {
                    // If there is already a Region table, load the region types that it has.
                    // We will reuse them if we write more regions of the same type, or add
                    // to the list if we get new RegionTypes
                    LoadRegionTypes();
                }

                // Register the standard region types
                if (!regionTypeCodes.ContainsKey(AdapterRegionType))
                    RegisterRegionType(AdapterRegionType);

                if (!regionTypeCodes.ContainsKey(InsertRegionType))
                    RegisterRegionType(InsertRegionType);

                if (!regionTypeCodes.ContainsKey(HqRegionType))
                    RegisterRegionType(HqRegionType);
            }

            /// <summary>
            /// If we are appending to an existing Region table, read in the existing region types
            /// </summary>
            private void LoadRegionTypes()
            {
                var typeNamesAttr = regionDataset.GetAttribute(TypesAttribute);

                if (typeNamesAttr == null)
                    return;

                var regionType = (string[]) regionDataset.GetAttribute(TypesAttribute).Read();
                var regionDescription = (string[]) regionDataset.GetAttribute(DescriptionAttribute).Read();
                var source = (string[]) regionDataset.GetAttribute(SourceAttribute).Read();

                for (int i = 0; i < regionType.Length; i++)
                {
                    var rt = new RegionType
                        {
                            Description = regionDescription[i],
                            Source = source[i],
                            Type = regionType[i]
                        };

                    var regionIdx = regionTypeList.Count;
                    regionTypeList.Add(rt);
                    regionTypeCodes[rt] = regionIdx;
                }
            }

            // The Region table contains three string array attribute that form a simple 'RegionType table'
            // The columns of the table are 'RegionTypes', 'RegionDescriptions', and 'RegionSources'.
            // This method add a new region type to the list of region types in this table
            private void RegisterRegionType(RegionType newRegionType)
            {
                //  To register a new RegionType, we rewrite the complete RegionType table.
                var regionId = regionTypeList.Count();
                regionTypeList.Add(newRegionType);
                regionTypeCodes[newRegionType] = regionId;

                WriteReigionAttribute(TypesAttribute, r => r.Type);
                WriteReigionAttribute(DescriptionAttribute, r => r.Description);
                WriteReigionAttribute(SourceAttribute, r => r.Source);
            }

            /// <summary>
            /// Write attributes of the Region table containing on column of the RegionType table
            /// </summary>
            /// <param name="attributeName"></param>
            /// <param name="f"></param>
            private void WriteReigionAttribute(string attributeName, Func<RegionType, string> f)
            {
                var attrOut = regionTypeList.Map(f);
                regionDataset.InsertAttribute(attributeName, attrOut);
            }

            public void Write(IList<Region> regions)
            {
                var n = regions.Count;

                var dims = regionDataset.Dataspace.Dimensions;
                var startRow = dims[0];
                var width = dims[1];

                regionDataset.Extend(new[] {startRow + n, width});
                var newDs = regionDataset.Dataspace;

                newDs.SelectHyperslab(new[] {startRow, 0}, null, new[] {n, width}, null);
                var output = regions.Map(MakeRow).JaggedToRectRows();
                regionDataset.Write(output, newDs);
            }

            // Convert a Region structure to a row of the region table
            private int[] MakeRow(Region r)
            {
                if (r.Type == null)
                    throw new ArgumentException("No RegionType specified in Region spec");

                if (!regionTypeCodes.ContainsKey(r.Type))
                    RegisterRegionType(r.Type);


                return new[] {r.HoleNumber, regionTypeCodes[r.Type], r.Start, r.End, r.Score};
            }
        }
    }
}
