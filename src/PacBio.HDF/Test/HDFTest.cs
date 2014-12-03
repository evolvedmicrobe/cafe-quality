using System;
using System.Linq;
using System.Diagnostics;
using System.IO;
using System.Text;
using NUnit.Framework;

namespace PacBio.HDF.Test
{
    /// <summary>
    /// Unit tests for my HDF chunk stuff
    /// </summary>
    [TestFixture]
    public class HDFTest
    {

        // FIXME - use a temp filename
        private string filename = "high.h5";

        [Test]
        public void TestAttributeWrite()
        {
            var high = HDFFile.Open(filename, FileMode.Create, FileAccess.ReadWrite);
            int[] expected = new int[] {1, 2, 3};
            var group = high.CreateGroup("/Root");

            group.InsertAttribute("3Ints", expected);

            var actual = group.GetAttribute("3Ints").Read() as int[];
            
            Assert.AreEqual(expected, actual);

            high.Dispose();
        }

        /// <summary>
        /// Sample of how to use from matlab
        /// </summary>
        [Test]
        public void TestHighLevel()
        {
            IHighLevelChunks high = new HighLevelChunks(filename, true);

            string [] strings = { "dog", "cat", "fish" };
   
            high.WriteDataset("/Subgroup/AString", "kevin was here");
            high.WriteDataset("/AnInt", 5);
            high.WriteDataset("/Subgroup/StringArray", strings);
            high.WriteDataset("Subgroup/Array", new double[] { 1.1, 1.2, 1.3 });

            CreateHighLevelFrames(high);

            string[] strs = (string []) high.ReadDataset("Subgroup/StringArray");
            Console.WriteLine("readback " + strs[1]);
            Assert.AreEqual(strings[1], strs[1]);
            Console.WriteLine("readback " + high.ReadDataset("Subgroup/AString"));

            // Make sure missing things come back as null
            var res = high.ReadDataset("Fish/Dog");
            Assert.IsNull(res);
        }

        /// <summary>
        /// Test the hyperslab high level API
        /// </summary>
        /// <param name="high"></param>
        void CreateHighLevelFrames(IHighLevelChunks high)
        {
            // True if we want to only add frames as needed
            bool growDynamically = true;

            // For the high level API we always need to start with at least one entry filled in the dataset (so we can
            // dynamically detect matlab types)
            long[] dims = { growDynamically ? 1 : numFrames, height, width };
            long[] maxDims = { growDynamically ? -1 : numFrames, height, width };

            ushort[] frame = new ushort[height * width];

            // The frame we are writing now
            long[] start = { 0, 0, 0 };
            long[] stride = { 1, 1, 1 };          // Step by 1 in each dir
            long[] count = { 1, 1, 1 };           // We are only writing one block at a time
            long[] block = { 1, height, width };  // One block per frame

            // IDataspace dspace = high.File.CreateDataspace(
            //    dims,
            //    maxDims);

            // Put a nice pattern into our frame
            for (int fi = 0; fi < frame.Length; fi++)
                frame[fi] = (ushort)fi;

#pragma warning disable 219
            IDataset dset = null;
            for (int i = 0; i < numFrames; i++)
            {
#if false
                // no longer needed - we now auto extend on write
                if (growDynamically && dset != null)
                {
                    dims[0] = i + 1; // Update # frames
                    dset.Extend(dims);
                }
#endif
                IDataspace dspace = high.File.CreateDataspace(
                    dims,
                    maxDims);

                frame[0] = (ushort)i;

                // Set the frame we are writing to
                start[0] = i;
                dspace.SelectHyperslab(start, stride, count, block);

                dset = high.WriteDataset("Movie/Frames", frame, dspace);

                // It is important that we manipulate the dataspace associated with the dataset - otherwise hyperslabs
                // won't work
                // dspace = dset.Dataspace;
            }
        }

        static string GetSpaces(int numspace)
        {
            StringBuilder b = new StringBuilder();

            while (numspace-- >= 0)
                b.Append(" ");

            return b.ToString();
        }


        void ReadDataset(IChunkElement arg, int depth)
        {
            int numspace = depth * 2;

            Console.WriteLine(GetSpaces(numspace) + arg);

            var dset = arg as IDataset;
            if (dset != null)
            {
                object o = dset.Read();
                Console.WriteLine(GetSpaces(numspace) + "Data: " + o);
            }

            var targ = arg as IAttributeTarget;
            if (null != targ)
                foreach (IDataContainer a in targ.GetAttributes())
                {
                    Console.WriteLine(GetSpaces(numspace) + "  Attr: " + a + " = " + a.Read());
                    a.Dispose();        // No leaking
                }

            var group = arg as IGroup;
            if (group != null)
            {
                foreach (IChunkElement child in group.GetChildren())
                {
                    ReadDataset(child, depth + 1);
                    child.Dispose();
                }
            }
        }

        // FIXME - use a temp file
        private readonly string testFile = "test.h5";

        [Test]
        public void ReadFile()
        {
            using (IChunkFile file = HighLevelChunks.Open(testFile, FileMode.Open, FileAccess.Read))
            {
                using (IGroup group = (IGroup)file.GetChild("Movie"))
                {
                    ReadFrames(group);
                }

                ReadDataset(file, 0);
            }
        }


        /// <summary>
        /// Read frame subchunks from the file
        /// </summary>
        /// <param name="movie"></param>
        void ReadFrames(IGroup movie)
        {
            using (IDataset dset = (IDataset)movie.GetChild("Frames"))
            {
                ushort[] frame = new ushort[width * height];

                long[] start = { 0, 0, 0 };
                long[] stride = { 1, 1, 1 };          // Step by 1 in each dir
                long[] count = { 1, 1, 1 };           // We are only writing one block at a time
                long[] block = { 1, height, width };  // One block per frame

                Array array = frame;

                for (int i = 0; i < numFrames; i++)
                {
                    start[0] = i;

                    // FIXME - change Dataspace to not leak objects
                    using (IDataspace dspace = dset.Dataspace)
                    {
                        dspace.SelectHyperslab(start, stride, count, block);

                        dset.Read(ref array, dspace);
                    }

                    if (frame[0] != i)
                        Debug.Fail("frame read test failed on frame " + i);
                }

                
            }

            using (IDataset refSet = (IDataset)movie.GetChild("RefPtr"))
            {
                IReference refptr = (IReference) refSet.Read();
                Console.WriteLine("Got dset: " + refptr.Dereference());
                Console.WriteLine("Got dspace: " + refptr.GetRegion());
            }

            Console.WriteLine("Frame read test succeeds");
        }


        int width = 512;
        int height = 100;
        int numFrames = 10;

        /// <summary>
        /// Pretend to write the frames of a movie
        /// </summary>
        /// This test has shown that writing frames is about 50 frames/sec for astro sized frames (I don't remember what RIFF
        /// files were for comparison) when chunking/dynamic growing is on.
        /// When chunking is _off_ things are muck faster - about 200 frames/sec.
        private void CreateFrames(IGroup group)
        {
            // True if we want to only add frames as needed
            bool growDynamically = false;

            long[] dims = {growDynamically ? 0 : numFrames, height, width};
            long[] maxDims = {growDynamically ? -1 : numFrames, height, width};

            ushort[] frame = new ushort[height*width];

            // The frame we are writing now
            long[] start = {0, 0, 0};
            long[] stride = {1, 1, 1}; // Step by 1 in each dir
            long[] count = {1, 1, 1}; // We are only writing one block at a time
            long[] block = {1, height, width}; // One block per frame
            IReference refptr = null;

            using (IDataspace initialDspace = group.File.CreateDataspace(
                dims,
                maxDims))
            {
                // Put a nice pattern into our frame
                for (int fi = 0; fi < frame.Length; fi++)
                    frame[fi] = (ushort) fi;

                using (var ushortType = group.File.CreateDatatype(typeof (ushort)))
                {
                    using (IDataset dset = group.CreateDataset("Frames",
                                                               ushortType,
                                                               initialDspace))
                    {
                        for (int i = 0; i < numFrames; i++)
                        {
                            if (growDynamically)
                            {
                                dims[0] = i + 1; // Update # frames
                                dset.Extend(dims);
                            }

                            frame[0] = (ushort) i;

                            // Set the frame we are writing to
                            start[0] = i;

                            // FIXME - change Dataspace to not leak objects
                            using (IDataspace dspace = dset.Dataspace)
                            {
                                dspace.SelectHyperslab(start, stride, count, block);

                                // write our frame
                                dset.Write(frame, dspace);

                                if (refptr == null)
                                    // Create a test reference
                                    refptr = dset.CreateReference(dspace);
                            }
                        }

                    }
                }
            }

            using (var simplespace = group.File.CreateDataspace())
            {
                using (var dtype = group.File.CreateDatatype(refptr.GetType()))
                using (IDataset refSet = group.CreateDataset("RefPtr", dtype, simplespace))
                {
                    refSet.Write(refptr);
                }
            }
        }

        [Test]
        public void CreateFile()
        {
            using (IChunkFile file = HighLevelChunks.Open(testFile, FileMode.Create, FileAccess.ReadWrite))
            {
                using (IGroup group = file.CreateGroup("Movie"))
                {
                    CreateFrames(group);
                }


                using (var dspace = file.CreateDataspace())
                {
                    using (IDatatype inttype = file.CreateDatatype(typeof (int)))
                    {
                        using (IDataset dset = file.CreateDataset("simple", inttype, dspace))
                        {
                            using (
                                IDataContainer attr =
                                    dset.CreateAttribute("simpattr", file.CreateDatatype(typeof (string)), dspace))
                            {
                                attr.Write("I hate monkies");
                            }
                            dset.Write(5);
                        }
                    }

                    ChunkUtils.Create(file, "simpleCreate", "some string");
                    ChunkUtils.Create(file, "hierarchical/create1", "some other string 1");
                    ChunkUtils.Create(file, "hierarchical/create2", "some other string 2");

                    using (IGroup group = file.CreateGroup("group"))
                    {
                        using (IDatatype dtype = file.CreateDatatype(typeof (double)))
                        {
                            using (IDataset dset = group.CreateDataset("adouble", dtype, dspace))
                                dset.Write(4.20);
                        }

                        ChunkUtils.Create(file, "inagroup/create", 3.14);

                        using (IDatatype strtype = file.CreateDatatype(typeof (string)))
                        {
                            using (IDataset dset = group.CreateDataset("astr", strtype, dspace))
                                dset.Write("Kevin's string");
                        }

                        int[,] iarray = new int[,]
                            {
                                {1, 2},
                                {3, 4},
                                {5, 6}
                            };

                        using (IDataset dset = file.CreateDataset("iarray",
                                                                  file.CreateDatatype(typeof (int)),
                                                                  file.CreateDataspace(new long[] {3, 2},
                                                                                       new long[] {3, 2}))
                            )
                            dset.Write(iarray);

                        ushort[,] uarray = new ushort[,]
                            {
                                {1, 2},
                                {3, 4},
                                {5, 6}
                            };

                        using (IDataset dset = file.CreateDataset("uarray",
                                                                  file.CreateDatatype(typeof (ushort)),
                                                                  file.CreateDataspace(new long[] {3, 2},
                                                                                       new long[] {3, 2}))
                            )
                            dset.Write(uarray);

                        double[,] darray = new double[,]
                            {
                                {1.1, 1.2},
                                {2.1, 2.2},
                                {3.1, 3.2}
                            };

                        double[,] darray2 = new double[,]
                            {
                                {1.1, 1.2},
                                {2.1, 2.2},
                                {3.1, 3.2},
                                {4.1, 4.2}
                            };


                        using (IDataset dset = file.CreateDataset("array",
                                                                  file.CreateDatatype(typeof (double)),
                                                                  file.CreateDataspace(new long[] {3, 2},
                                                                                       new long[] {6, 6})))
                        {
                            dset.Write(darray);
                            dset.Extend(new long[] {4, 2});
                            dset.Write(darray2);
                        }
                    }
                }
            }
        }
    }
}
