using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace SuiteSparseCSWrap
{
    /// <summary>
    /// Sparse Matrix class. Defined using triplets format. 
    /// Can be square or non-square, upper/lower triangle or non-symmetric 
    /// </summary>
    public class SparseMatrix
    {
        private Dictionary< Tuple<int, int>, double> data = new Dictionary<Tuple<int,int>,double>();
        private int sym_type = 0; // lower_triangle(-1), upper_triangle(1), unsymmetric(0);
        private void Swap(ref int i, ref  int j){
            int k = i; i = j; j = k;
        }
        private int Nrows = 0, Ncols = 0;
        
        /// <summary>
        /// Sets a value to (i,j) matrix entry
        /// </summary>
        /// <param name="i">row index</param>
        /// <param name="j">column index</param>
        /// <param name="val">value</param>
        public void SetValue(int i, int j, double val)
        {
            add_dim(i, j);
            if (sym_type < 0 && i < j) Swap(ref i, ref j); 
            data[Tup(i, j)] = val;
        }
        /// <summary>
        /// Gets a value from (i, j)
        /// </summary>
        /// <param name="i">row index</param>
        /// <param name="j">column index</param>
        /// <returns>value</returns>
        public double GetValue(int i, int j)
        {
            if (sym_type < 0 && i < j) Swap(ref i, ref j); 
            var ind = Tup(i, j);
            return data.ContainsKey(ind) ? data[ind] : 0.0;
        }
        /// <summary>
        /// returns all values as a list. Only Lower/upper parts for symmetrical matricies
        /// </summary>
        /// <returns>list of values (row index, column index, value) </returns>
        public List<Tuple<int, int, double>> GetAllValues()
        {
            List<Tuple<int, int, double>> ret = new List<Tuple<int, int, double>>();
            foreach (var e in data)
            {
                ret.Add(new Tuple<int, int, double>(e.Key.Item1, e.Key.Item2, e.Value));
            }
            return ret;
        }

        /// <summary>
        /// Clears all values
        /// </summary>
        public void Clear()
        {
            data.Clear();
            Ncols = Nrows = 0;
        }

        /// <summary>
        /// Get number of defined entries
        /// </summary>
        /// <returns>Number of entries</returns>
        public int GetNnz()
        {
            return data.Count;
        }
        
        
        /// <summary>
        /// Matrix dimension
        /// </summary>
        /// <returns> (number of rows, number of columns) </returns>
        public Tuple<int, int> GetDim()
        {
            return new Tuple<int, int>(Nrows, Ncols);
        }
        private void add_dim(int i, int j)
        {
            if (sym_type != NONSYM)
            {
                i = j = Math.Max(i, j);
            }
            if (i + 1 > Nrows) Nrows = i + 1;
            if (j + 1 > Ncols) Ncols = j + 1;
        }
        private void set_dimensions()
        {
            Ncols = Nrows = 0;
            foreach(var d in data) add_dim(d.Key.Item1, d.Key.Item2);
        }

        /// <summary>
        /// Lower triangle indicator
        /// </summary>
        public const int LOWERTRI = -1;
        /// <summary>
        /// Upper Triangle indicator
        /// </summary>
        public const int UPPERTRI = 1;
        /// <summary>
        /// Nonsym matrix indicator
        /// </summary>
        public const int NONSYM = 0;
        /// <summary>
        /// Symmetric type
        /// </summary>
        /// <returns>Possible values: LOWERTRI, UPPERTRI, NONSYM</returns>
        public int GetSymType()
        {
            return sym_type;
        }
        /// <summary>
        /// Set Symmetric type of Matrix (Lower Triangle, Upper Triangle, NonSymmetric)
        /// if Nonsymmetric matrix was changed to Triangle, all values outside defined triangles will be lost.
        /// if Triangle matrix was change to Nonsymmetric its (i,j) values will be copied to (j, i).
        /// So the Number of entries will be doubled.
        /// </summary>
        /// <param name="tp">Possible values: LOWERTRI, UPPERTRI, NONSYM</param>
        public void SetSymmetricType(int tp)
        {
            if (tp == sym_type) return;
            Dictionary<Tuple<int, int>, double> data2 = new Dictionary<Tuple<int, int>, double>();
            if (sym_type == 0)
            {
                if (tp < 0)  //unsym -> lower
                {
                    foreach (var e in data)
                        if (e.Key.Item1 >= e.Key.Item2) data2.Add(e.Key, e.Value);
                }
                else if (tp > 0) //unsym -> upper
                {
                    foreach (var e in data)
                        if (e.Key.Item1 <= e.Key.Item2) data2.Add(e.Key, e.Value);
                }
            } else if (sym_type == 1) {
                if (tp < 0) //upper -> lower
                {
                    foreach (var e in data) data2.Add(new Tuple<int, int>(e.Key.Item2, e.Key.Item1), e.Value);
                }
                else if (tp == 0) //upper -> unsym
                {
                    foreach (var e in data)
                    {
                        data2.Add(e.Key, e.Value);
                        data2.Add(new Tuple<int, int>(e.Key.Item2, e.Key.Item1), e.Value);
                    }
                }
            }
            else if (sym_type == -1) 
            {
                if (tp > 0) //lower->upper
                {
                    foreach (var e in data) data2.Add(new Tuple<int, int>(e.Key.Item2, e.Key.Item1), e.Value);
                }
                else if (tp == 0) //lower -> nonsym
                {
                    foreach (var e in data)
                    {
                        data2.Add(e.Key, e.Value);
                        data2.Add(new Tuple<int, int>(e.Key.Item2, e.Key.Item1), e.Value);
                    }
                }
            }
            sym_type = tp;
            data = data2;
            set_dimensions();
        }

        /// <summary>
        /// Matrix-Vector multiplication
        /// </summary>
        /// <param name="x">vector</param>
        /// <returns>Mat*x</returns>
        public double[] MultVect(double[] x)
        {
            double[] ret = new double[Nrows];
            foreach (var e in data)
            {
                ret[e.Key.Item1] += x[e.Key.Item2]*e.Value;
                if (sym_type != 0 && e.Key.Item1 != e.Key.Item2)
                {
                    ret[e.Key.Item2] += x[e.Key.Item1] * e.Value;
                }
            }
            return ret;
        }

        /// <summary>
        /// (Mat*x - rhs) maximum norm
        /// </summary>
        /// <param name="rhs">right hand side of the linear system</param>
        /// <param name="x">solution of the linear system</param>
        /// <returns>||Mat*x-rhs||</returns>
        public double MaxNorm(double[] rhs, double[] x)
        {
            double[] Ax = MultVect(x);
            double ret = 0;
            for (int i = 0; i < rhs.Length; ++i)
            {
                double v = Math.Abs(rhs[i] - Ax[i]);
                if (v > ret) ret = v;
            }
            return ret;
        }

        /// <summary>
        /// WriteMatrix into a MatrixMarket format
        /// </summary>
        /// <param name="fn">filename</param>
        public void WriteToFile(string fn)
        {
            //header
            System.IO.StreamWriter file = new System.IO.StreamWriter(fn);
            string tp = (sym_type == 0) ? "general" : "symmetric";
            file.WriteLine("%%MatrixMarket matrix coordinate real {0}", tp);
            var dim = GetDim();
            var Nnz = GetNnz();
            file.WriteLine("{0} {1} {2}", dim.Item1, dim.Item2, Nnz); 
            //elements
            foreach (var e in data)
            {
                if (sym_type == 0){
                    file.WriteLine("{0} {1} {2}", e.Key.Item1 + 1, e.Key.Item2 + 1, e.Value);
                }
                else if (sym_type < 0)
                {
                    if (e.Key.Item1 >= e.Key.Item2)
                    {
                        file.WriteLine("{0} {1} {2}", e.Key.Item1 + 1, e.Key.Item2 + 1, e.Value);
                    }
                }
                else if (sym_type > 0)
                {
                    if (e.Key.Item1 <= e.Key.Item2)
                    {
                        file.WriteLine("{0} {1} {2}", e.Key.Item2 + 1, e.Key.Item1 + 1, e.Value);
                    }
                }
                
            }
            file.Close();
        }
        
        /// <summary>
        /// Fill from MatrixMarket format file. All presenting entries will be freed.
        /// </summary>
        /// <param name="fn">filename</param>
        public void FillFromFile(string fn)
        {
            Clear();
            string[] lines = System.IO.File.ReadAllLines(fn);
            //1)header
            int i = 0;
            while (i>=lines.Length || lines[i].Trim().Length < 2 || lines[i].Trim()[0] != '%' || lines[i].Trim()[1] != '%') ++i;
            if (i == lines.Length)
                throw new Exception("Cannot read MatrixMarket file");
            
            string[] header = lines[i].Split();
            if (header.Length < 5 || header[0] != "%%MatrixMarket" || header[1] != "matrix" || header[2] != "coordinate" || header[3] != "real")
                throw new Exception("Cannot read MatrixMarket file");
            if (header[4] == "symmetric") sym_type = -1;
            else sym_type = 0;

            //2) m, n values
            ++i;
            int M=0, N=0, Nnz=0;
            while (i < lines.Length)
            {
                string line = lines[i].Trim();
                if (line.Length != 0 && line[0] != '%')
                {
                    string[] nums = line.Split();
                    if (nums.Length >= 3)
                    {
                        M = int.Parse(nums[0]);
                        N = int.Parse(nums[1]);
                        Nnz = int.Parse(nums[2]);
                        break;
                    }
                }
                ++i;
                if (i == lines.Length)
                    throw new Exception("Cannot read MatrixMarket file");
            }
            //fill values
            ++i;
            while (i < lines.Length || GetNnz() < Nnz)
            {
                string line = lines[i].Trim();
                if (line.Length != 0 && line[0] != '%')
                {
                    string[] nums = line.Split();
                    if (nums.Length >= 3)
                    {
                        int Row = int.Parse(nums[0]);
                        int Col = int.Parse(nums[1]);
                        double Val = double.Parse(nums[2]);
                        SetValue(Row-1, Col-1, Val);
                    }
                }
                ++i;
            }
            if (GetNnz()!=Nnz)
                throw new Exception("Cannot read MatrixMarket file");
        }

        private Tuple<int, int> Tup(int i, int j)
        {
            return new Tuple<int, int>(i, j);
        }
    }
}
