using System;
using System.IO;
using System.Linq;
using NUnit.Framework;
using PacBio.Utils;

namespace PacBio.HDF.Test
{
    public class TestType
    {
        public int a
        {
            get; set;
        }

        public float b
        {
            get;
            set;
        }


        public double[] c
        {
            get;
            set;
        }

        public int[] d
        {
            get;
            set;
        }


        public static TestType[] TestValues()
        {
            var v1 = new TestType { a = 5, b = -1.0f, c = new[] { 1.0, 2, 3 }, d = new[] { 1, 2, 3 } };
            var v2 = new TestType { a = -1, b = 1.37f, c = new[] { 1.0, 2, 3 }, d = new[] { 1, 2, 3 } };
            var v3 = new TestType { a = 1231, b = 6.0f, c = new[] { 1.0, 2, 3, 5, 6, 7, 8 }, d = new[] { 1, 2, 3, 5, 6, 7, 8 } };
            var v4 = new TestType { a = -123,   b = 6.0f, c = new double[0], d = new int[0] };
            var v5 = new TestType { a = 8,      b = 6.0f, c = new[] { 1.0 }, d = new[] {5}};

            return new[] {v1, v2, v3, v4, v5};
        }
    }


    public class TestTypeWriter : ArrayRecordStore<TestType>.Writer
    {
        public TestTypeWriter(IGroup group) : base(group)
        {
            WriteSetupSingleton("a", t => t.a);
            WriteSetupSingleton("b", t => t.b);

            var arrGroup = MakeParallelArrayGroup("Index");
            arrGroup.AddArrayField("c", t => t.c);
            arrGroup.AddArrayField("d", t => t.d);
            arrGroup.Close();
        }
    }


    public class TestTypeReader : ArrayRecordStore<TestType>.Reader
    {
        private Func<int, TestType> reader;

        public TestTypeReader(IGroup group)
            : base(group)
        {
            var ar = MakeSingletonReader<int>("a");
            var br = MakeSingletonReader<float>("b");

            var cr = MakeArrayReader<double>("c");
            var dr = MakeArrayReader<int>("d");

            reader =
                n => {
                    return new TestType
                               {
                                   a = ar.Read(n, 1)[0],
                                   b = br.Read(n, 1)[0],
                                   c = cr.Read(n, 1)[0],
                                   d = dr.Read(n, 1)[0]
                               };
                };
        }

        public TestType[] Read()
        {
            return NumRecords.Fill(i => reader(i));
        }

    }


    [TestFixture,Explicit]
    public class ArrayRecordStore
    {
        string testUri = "ArrayRecordStoreTest.h5";

        [Test]
        public void WriteData()
        {
            using (var hlc = HDFFile.Open(testUri, FileMode.Open, FileAccess.ReadWrite))
            {
                var group = hlc.CreateGroup("Test");
                var writer = new TestTypeWriter(group);

                writer.WriteRecords(TestType.TestValues());
            }
        }

        [Test]
        public void ReadData()
        {
            WriteData();

            // FIXME 
            /* 

            using (var hlc = new HighLevelChunks(testUri, true))
            {
                var group = hlc.CreateGroup("Test");
                var reader = new TestTypeReader(group);

                var answers = TestType.TestValues();
                var data = reader.Read();
                

                Assert.AreEqual(answers.Length, data.Length);

                for (int i = 0; i < answers.Length; i++)
                {
                    Assert.AreEqual(answers[i].a, answers[i].a);
                    Assert.AreEqual(answers[i].b, answers[i].b);
                    Assert.That(answers[i].c.SequenceEqual(data[i].c));
                    Assert.That(answers[i].d.SequenceEqual(data[i].d));
                }
            }
             */
        }
    }
}
