using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using PacBio.Utils;

namespace PacBio.IO
{


    public class CsvWriter : IDisposable
    {
        public CsvWriter(string filename)
        {
            f = new FileStream(filename, FileMode.Create);
            w = new StreamWriter(f);

            curRow = new object[0];
            colNames = new string[0];
        }

        private readonly FileStream f;
        private readonly StreamWriter w;
        private string[] colNames;
        private object[] curRow;
        private bool firstRow = true;
        private readonly Dictionary<string, int> columnIndicies = new Dictionary<string, int>();


        public static T[] Append<T>(T[] arr, T item)
        {
            var arrNew = new T[arr.Length + 1];
            arr.CopyTo(arrNew, 0);
            arrNew[arr.Length] = item;
            return arrNew;
        }

        private void addCol(string name, object data)
        {
            var n = curRow.Length;
            columnIndicies[name] = n;
            curRow = Append(curRow,data);
            colNames = Append(colNames, name);
        }

        public void AddCol(string name, object data)
        {
            int col;
            bool colExists = columnIndicies.TryGetValue(name, out col);

            if (colExists)
                curRow[col] = data;
            else
                addCol(name, data);
        }

        public object this[string col]
        {
            set { AddCol(col, value); }
        }

        public void Row()
        {
            if (firstRow)
            {
                csvLine(colNames.Cast<object>().ToArray());
                firstRow = false;
            }

            csvLine(curRow);
        }

        private void csvLine(object[] vals)
        {
            w.WriteLine(String.Join(",", vals.Map(v => v.ToString())));
        }

        private bool disposed = false;

        protected virtual void Dispose(bool disposing)
        {
            if (disposed)
                return;

            w.Close();
            f.Close();

            if (disposing)
            {
                w.Dispose();
                f.Dispose();
            }

            disposed = true;
        }

        public void Dispose()
        {
            Dispose(true);
            GC.SuppressFinalize(this);
        }

        ~CsvWriter()
        {
            Dispose(false);
        }
    }
}