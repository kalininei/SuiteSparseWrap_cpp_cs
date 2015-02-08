using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using SuiteSparseCSWrap;

namespace SolverTest
{
    class Program
    {

        // =============== Random number generators:
        static System.Random r = new System.Random(0);
        static int get_random_int(int min, int max)
        {
            return r.Next(min, max + 1);
        }
        static double get_random_double(Tuple<double, double> bnd)
        {
            return r.NextDouble() * (bnd.Item2 - bnd.Item1) + bnd.Item1;
        }
        
        // ================= Tests
        //diagonal matrix solver
        static void test1()
        {
            int N = 100;
            SparseMatrix A = new SparseMatrix();
            double[] rhs = new double[N];
            double[] x = new double[N];
            for (int i = 0; i < N; ++i)
            {
                A.SetValue(i, i, 1.0);
                rhs[i] = 12 * i + Math.Sin(i);
            }
            Solver Slv = new Solver(A, Solver.CHOLMODSOLVER);
            Slv.Solve(rhs, x);
            //X==RHS
            double diff = 0;
            for (int i = 0; i < N; ++i)
            {
                if (Math.Abs(rhs[i] - x[i]) > diff) diff = Math.Abs(rhs[i] - x[i]);
            }

            if (Math.Abs(diff) < 1e-12)
            {
                Console.WriteLine("Test: passed");
            }
            else
            {
                Console.WriteLine("Test: failed");
            }

            //IMplicit memory clean. Could be omitted.
            Slv.Dispose();
        }

        //square nonsymmetric matrix solver cholmod solver
        static void test2()
        {
            SparseMatrix A = new SparseMatrix();
            double[] rhs = {1, 0, 0, 1};
            double[] x = new double[4];
            A.SetValue(0, 0, 2.0); A.SetValue(0, 3, 4.0);
            A.SetValue(1, 2, 1.0); A.SetValue(1, 3, 2.0);
            A.SetValue(2, 0, 12); A.SetValue(2, 1, 3.0); A.SetValue(2, 2, 7.0);
            A.SetValue(3, 2, 2.0); A.SetValue(3, 3, 1.0);
            Solver S = new Solver(A, Solver.CHOLMODSOLVER);
            S.Solve(rhs, x);
            //Norm
            double nrm = A.MaxNorm(rhs, x);
            if (Math.Abs(nrm) < 1e-12)
            {
                Console.WriteLine("Test: passed");
            }
            else
            {
                Console.WriteLine("Test: failed");
            }
        }

        //square nonsymmetric matrix solver umfpack solver
        static void test3()
        {
            SparseMatrix A = new SparseMatrix();
            double[] rhs = { 1, 0, 0, 1 };
            double[] x = new double[4];
            A.SetValue(0, 0, 2.0); A.SetValue(0, 3, 4.0);
            A.SetValue(1, 2, 1.0); A.SetValue(1, 3, 2.0);
            A.SetValue(2, 0, 12); A.SetValue(2, 1, 3.0); A.SetValue(2, 2, 7.0);
            A.SetValue(3, 2, 2.0); A.SetValue(3, 3, 1.0);
            Solver S = new Solver(A, Solver.UMFPACKSOLVER);
            S.Solve(rhs, x);
            //Norm
            double nrm = A.MaxNorm(rhs, x);
            if (Math.Abs(nrm) < 1e-12)
            {
                Console.WriteLine("Test: passed");
            }
            else
            {
                Console.WriteLine("Test: failed");
            }
        }

        //square nonsymmetric matrix qr factorization solver
        static void test4()
        {
            SparseMatrix A = new SparseMatrix();
            double[] rhs = { 1, 0, 0, 1 };
            double[] x = new double[4];
            A.SetValue(0, 0, 2.0); A.SetValue(0, 3, 4.0);
            A.SetValue(1, 2, 1.0); A.SetValue(1, 3, 2.0);
            A.SetValue(2, 0, 12); A.SetValue(2, 1, 3.0); A.SetValue(2, 2, 7.0);
            A.SetValue(3, 2, 2.0); A.SetValue(3, 3, 1.0);
            Solver S = new Solver(A, Solver.SPQRSOLVER);
            S.Solve(rhs, x);
            //Norm
            double nrm = A.MaxNorm(rhs, x);
            if (Math.Abs(nrm) < 1e-12)
            {
                Console.WriteLine("Test: passed");
            }
            else
            {
                Console.WriteLine("Test: failed");
            }
        }

        //non square matrix qr factorization solver
        static void test5()
        {
            // System is  x+y+z=1;
            //            x=0.5;
            //            y=0.5;
            //            z=1e-3

            SparseMatrix A = new SparseMatrix();
            double[] rhs = { 1, 0.5, 0.5, 1e-3 };
            double[] x = new double[4];
            A.SetValue(0, 0, 1.0); A.SetValue(0, 1, 1.0); A.SetValue(0, 2, 1.0);
            A.SetValue(1, 0, 1.0);
            A.SetValue(2, 1, 1.0);
            A.SetValue(3, 2, 1.0);
            Solver S = new Solver(A, Solver.SPQRSOLVER);
            S.Solve(rhs, x);
            //Norm
            double nrm = A.MaxNorm(rhs, x);
            if (Math.Abs(nrm) < 1e-3)
            {
                Console.WriteLine("Test: passed");
            }
            else
            {
                Console.WriteLine("Test: failed");
            }
        }

        //symmetric matrix tests
        static void test6()
        {
            SparseMatrix Mat = new SparseMatrix();
            Mat.SetSymmetricType(SparseMatrix.UPPERTRI);
            Mat.SetValue(0, 0, 2); Mat.SetValue(0, 3, 1);
            Mat.SetValue(1, 1, 1); Mat.SetValue(2, 2, 1);
            Mat.SetValue(3, 3, 5);
            double[] b = { 1, 1, 1, 1 };
            double[] x = new double[4];
            //CHOLMOD
            Solver CholS = new Solver(Mat, Solver.CHOLMODSOLVER);
            CholS.Solve(b, x);
            double CholNrm = Mat.MaxNorm(b, x);
            //UMFPACK
            Solver UMFS = new Solver(Mat, Solver.UMFPACKSOLVER);
            UMFS.Solve(b, x);
            double UMFNrm = Mat.MaxNorm(b, x);
            //QR
            Solver QRS = new Solver(Mat, Solver.SPQRSOLVER);
            QRS.Solve(b, x);
            double QRNrm = Mat.MaxNorm(b, x);
            if (QRNrm < 1e-12 && UMFNrm < 1e-12 && CholNrm < 1e-12)
            {
                Console.WriteLine("Test: passed");
            }
            else
            {
                Console.WriteLine("Test: failed");
            }

        }

        //upper trianlge positive definite matrix with band stencil
        static void test7(int Nrow, int[] Stencil, Tuple<double, double> valBnd)
        {
            //1) build matrix
            SparseMatrix Mat = new SparseMatrix();
            Mat.SetSymmetricType(SparseMatrix.UPPERTRI);
            //random non-diagonal
            for (int i = 0; i < Nrow; ++i)
            {
                foreach(int j in Stencil){
                    int ind = i + j;
                    if (ind < Nrow)
                    {
                        Mat.SetValue(i, ind, -get_random_double(valBnd));
                    }
                }
            }
            //diagonal value >= Sum(row)
            double[] diag = new double[Nrow];
            for (int i = 0; i < diag.Length; ++i ) diag[i] = get_random_double(valBnd);
            foreach (var entry in Mat.GetAllValues())
            {
                diag[entry.Item1] -= entry.Item3;
                diag[entry.Item2] -= entry.Item3;
            }
            for (int i = 0; i < diag.Length; ++i)
            {
                Mat.SetValue(i, i, diag[i]);
            }
            //2) build rhs
            double[] rhs1 = new double[Nrow];
            double[] rhs2 = new double[Nrow];
            double[] x = new double[Nrow];
            for (int i = 0; i < Nrow; ++i)
            {
                rhs1[i] = get_random_double(valBnd);
                rhs2[i] = get_random_double(valBnd);
            }

            // ------------- Cholezky solver
            Solver CholSlv = new Solver(Mat, Solver.CHOLMODSOLVER);

            //cold start time
            var watch1 = System.Diagnostics.Stopwatch.StartNew();
            CholSlv.Solve(rhs1, x);
            watch1.Stop();
            var CholTime1 = watch1.ElapsedMilliseconds;

            double Chol_nrm1 = Mat.MaxNorm(rhs1, x);
            
            //warm start time
            var watch2 = System.Diagnostics.Stopwatch.StartNew();
            CholSlv.Solve(rhs2, x);
            watch2.Stop();
            var CholTime2 = watch2.ElapsedMilliseconds;

            double Chol_nrm2 = Mat.MaxNorm(rhs2, x);


            // --------------- UMFPACK solver
            Solver UMFSlv = new Solver(Mat, Solver.UMFPACKSOLVER);
            //cold start time
            var watch3 = System.Diagnostics.Stopwatch.StartNew();
            UMFSlv.Solve(rhs1, x);
            watch3.Stop();
            var UMFTime1 = watch3.ElapsedMilliseconds;
            double UMF_nrm1 = Mat.MaxNorm(rhs1, x);

            //warm start time
            var watch4 = System.Diagnostics.Stopwatch.StartNew();
            UMFSlv.Solve(rhs2, x);
            watch4.Stop();
            var UMFTime2 = watch4.ElapsedMilliseconds;
            double UMF_nrm2 = Mat.MaxNorm(rhs2, x);


            // --------------- QR solver
            Solver QRSlv = new Solver(Mat, Solver.SPQRSOLVER);
            //cold start time
            var watch5 = System.Diagnostics.Stopwatch.StartNew();
            QRSlv.Solve(rhs1, x);
            watch5.Stop();
            var QRTime1 = watch5.ElapsedMilliseconds;
            double QR_nrm1 = Mat.MaxNorm(rhs1, x);

            //warm start time
            var watch6 = System.Diagnostics.Stopwatch.StartNew();
            QRSlv.Solve(rhs2, x);
            watch6.Stop();
            var QRTime2 = watch6.ElapsedMilliseconds;
            double QR_nrm2 = Mat.MaxNorm(rhs2, x);

            Console.WriteLine("Test results:");
            Console.WriteLine("\tCholmod: Cold Start {0} ms, Warm Start {1} ms, Norm {2}", CholTime1, CholTime2, Chol_nrm1);
            Console.WriteLine("\tUMFPACK: Cold Start {0} ms, Warm Start {1} ms, Norm {2}", UMFTime1, UMFTime2, UMF_nrm1);
            Console.WriteLine("\tQRSOLVER: Cold Start {0} ms, Warm Start {1} ms, Norm {2}", QRTime1, QRTime2, QR_nrm1);
        }

        //matrix with random stencil
        static void test8(int Nrow, int RowEntries, Tuple<double, double> valBnd)
        {
            //1) build matrix
            SparseMatrix Mat = new SparseMatrix();
            for (int i = 0; i < Nrow; ++i)
            {
                //diagonal
                Mat.SetValue(i, i, get_random_double(valBnd));
                //random non-diagonal
                for (int j = 0; j < RowEntries - 1; ++j)
                {
                    Mat.SetValue(i, get_random_int(0, Nrow-1), get_random_double(valBnd));
                }
            }

            //2) build rhs
            double[] rhs1 = new double[Nrow];
            double[] rhs2 = new double[Nrow];
            double[] x = new double[Nrow];
            for (int i = 0; i < Nrow; ++i)
            {
                rhs1[i] = get_random_double(valBnd);
                rhs2[i] = get_random_double(valBnd);
            }

            // ------------- Cholezky solver
            Solver CholSlv = new Solver(Mat, Solver.CHOLMODSOLVER);
            
            //cold start time
            var watch1 = System.Diagnostics.Stopwatch.StartNew();
            CholSlv.Solve(rhs1, x);
            watch1.Stop();
            var CholTime1 = watch1.ElapsedMilliseconds;
            
            double Chol_nrm1 = Mat.MaxNorm(rhs1, x);

            //warm start time
            var watch2 = System.Diagnostics.Stopwatch.StartNew();
            CholSlv.Solve(rhs2, x);
            watch2.Stop();
            var CholTime2 = watch2.ElapsedMilliseconds;
            
            double Chol_nrm2 = Mat.MaxNorm(rhs2, x);

            
            // --------------- UMFPACK solver
            Solver UMFSlv = new Solver(Mat, Solver.UMFPACKSOLVER);
            //cold start time
            var watch3 = System.Diagnostics.Stopwatch.StartNew();
            UMFSlv.Solve(rhs1, x);
            watch3.Stop();
            var UMFTime1 = watch3.ElapsedMilliseconds;
            double UMF_nrm1 = Mat.MaxNorm(rhs1, x);

            //warm start time
            var watch4 = System.Diagnostics.Stopwatch.StartNew();
            UMFSlv.Solve(rhs2, x);
            watch4.Stop();
            var UMFTime2 = watch4.ElapsedMilliseconds;
            double UMF_nrm2 = Mat.MaxNorm(rhs2, x);

            // --------------- QR solver
            Solver QRSlv = new Solver(Mat, Solver.SPQRSOLVER);
            //cold start time
            var watch5 = System.Diagnostics.Stopwatch.StartNew();
            QRSlv.Solve(rhs1, x);
            watch5.Stop();
            var QRTime1 = watch5.ElapsedMilliseconds;
            double QR_nrm1 = Mat.MaxNorm(rhs1, x);

            //warm start time
            var watch6 = System.Diagnostics.Stopwatch.StartNew();
            QRSlv.Solve(rhs2, x);
            watch6.Stop();
            var QRTime2 = watch6.ElapsedMilliseconds;
            double QR_nrm2 = Mat.MaxNorm(rhs2, x);

            // --------------- Dense solver
            Solver LUSlv = new Solver(Mat, Solver.DENSESOLVER);
            //cold start time
            var watch7 = System.Diagnostics.Stopwatch.StartNew();
            LUSlv.Solve(rhs1, x);
            watch7.Stop();
            var LUTime1 = watch7.ElapsedMilliseconds;
            double LU_nrm1 = Mat.MaxNorm(rhs1, x);

            //warm start time
            var watch8 = System.Diagnostics.Stopwatch.StartNew();
            LUSlv.Solve(rhs2, x);
            watch8.Stop();
            var LUTime2 = watch8.ElapsedMilliseconds;
            double LU_nrm2 = Mat.MaxNorm(rhs2, x);


            Console.WriteLine("Test results:");
            Console.WriteLine("\tCholmod: Cold Start {0} ms, Warm Start {1} ms, Norm {2}", CholTime1, CholTime2, Chol_nrm1);
            Console.WriteLine("\tUMFPACK: Cold Start {0} ms, Warm Start {1} ms, Norm {2}", UMFTime1, UMFTime2, UMF_nrm1);
            Console.WriteLine("\tQRSOLVER: Cold Start {0} ms, Warm Start {1} ms, Norm {2}", QRTime1, QRTime2, QR_nrm1);
            Console.WriteLine("\tDENSESOLVER: Cold Start {0} ms, Warm Start {1} ms, Norm {2}", LUTime1, LUTime2, LU_nrm1);

        }

        //Unsymmetric matrix with band stencil 
        static void test9(int Nrow, int[] Stencil, Tuple<double, double> valBnd)
        {
            //1) build matrix
            SparseMatrix Mat = new SparseMatrix();
            //random non-diagonal
            for (int i = 0; i < Nrow; ++i)
            {
                foreach (int j in Stencil)
                {
                    int ind1 = i + j, ind2 = i - j;
                    if (ind1 < Nrow)
                    {
                        Mat.SetValue(i, ind1, -get_random_double(valBnd));
                    }
                    if (ind2 >= 0)
                    {
                        Mat.SetValue(i, ind2, -get_random_double(valBnd));
                    }
                }
            }
            //diagonal value >= Sum(row)
            double[] diag = new double[Nrow];
            for (int i = 0; i < diag.Length; ++i) diag[i] = get_random_double(valBnd);
            foreach (var entry in Mat.GetAllValues())
            {
                diag[entry.Item1] -= entry.Item3;
                diag[entry.Item2] -= entry.Item3;
            }
            for (int i = 0; i < diag.Length; ++i)
            {
                Mat.SetValue(i, i, diag[i]);
            }
            //2) build rhs
            double[] rhs1 = new double[Nrow];
            double[] rhs2 = new double[Nrow];
            double[] x = new double[Nrow];
            for (int i = 0; i < Nrow; ++i)
            {
                rhs1[i] = get_random_double(valBnd);
                rhs2[i] = get_random_double(valBnd);
            }

            // ------------- Cholezky solver
            Solver CholSlv = new Solver(Mat, Solver.CHOLMODSOLVER);

            //cold start time
            var watch1 = System.Diagnostics.Stopwatch.StartNew();
            CholSlv.Solve(rhs1, x);
            watch1.Stop();
            var CholTime1 = watch1.ElapsedMilliseconds;

            double Chol_nrm1 = Mat.MaxNorm(rhs1, x);

            //warm start time
            var watch2 = System.Diagnostics.Stopwatch.StartNew();
            CholSlv.Solve(rhs2, x);
            watch2.Stop();
            var CholTime2 = watch2.ElapsedMilliseconds;

            double Chol_nrm2 = Mat.MaxNorm(rhs2, x);


            // --------------- UMFPACK solver
            Solver UMFSlv = new Solver(Mat, Solver.UMFPACKSOLVER);
            //cold start time
            var watch3 = System.Diagnostics.Stopwatch.StartNew();
            UMFSlv.Solve(rhs1, x);
            watch3.Stop();
            var UMFTime1 = watch3.ElapsedMilliseconds;
            double UMF_nrm1 = Mat.MaxNorm(rhs1, x);

            //warm start time
            var watch4 = System.Diagnostics.Stopwatch.StartNew();
            UMFSlv.Solve(rhs2, x);
            watch4.Stop();
            var UMFTime2 = watch4.ElapsedMilliseconds;
            double UMF_nrm2 = Mat.MaxNorm(rhs2, x);


            // --------------- QR solver
            Solver QRSlv = new Solver(Mat, Solver.SPQRSOLVER);
            //cold start time
            var watch5 = System.Diagnostics.Stopwatch.StartNew();
            QRSlv.Solve(rhs1, x);
            watch5.Stop();
            var QRTime1 = watch5.ElapsedMilliseconds;
            double QR_nrm1 = Mat.MaxNorm(rhs1, x);

            //warm start time
            var watch6 = System.Diagnostics.Stopwatch.StartNew();
            QRSlv.Solve(rhs2, x);
            watch6.Stop();
            var QRTime2 = watch6.ElapsedMilliseconds;
            double QR_nrm2 = Mat.MaxNorm(rhs2, x);

            Console.WriteLine("Test results:");
            Console.WriteLine("\tCholmod: Cold Start {0} ms, Warm Start {1} ms, Norm {2}", CholTime1, CholTime2, Chol_nrm1);
            Console.WriteLine("\tUMFPACK: Cold Start {0} ms, Warm Start {1} ms, Norm {2}", UMFTime1, UMFTime2, UMF_nrm1);
            Console.WriteLine("\tQRSOLVER: Cold Start {0} ms, Warm Start {1} ms, Norm {2}", QRTime1, QRTime2, QR_nrm1);
        }

        //Unsymmetric matrix with band stencil 
        static void test10(int Nrow, Tuple<double, double> valBnd)
        {
            //suppress warning
            Solver.NoWarnings = true;
            //1) build matrix
            SparseMatrix Mat = new SparseMatrix();
            //fill with random values
            for (int i = 0; i < Nrow; ++i)
            {
                for (int j = 0; j < Nrow; ++j)
                {
                    Mat.SetValue(i, j,get_random_double(valBnd));
                }
            }

            //2) build rhs
            double[] rhs1 = new double[Nrow];
            double[] rhs2 = new double[Nrow];
            double[] x = new double[Nrow];
            for (int i = 0; i < Nrow; ++i)
            {
                rhs1[i] = get_random_double(valBnd);
                rhs2[i] = get_random_double(valBnd);
            }

            // ------------- Cholezky solver
            Solver CholSlv = new Solver(Mat, Solver.CHOLMODSOLVER);

            //cold start time
            var watch1 = System.Diagnostics.Stopwatch.StartNew();
            CholSlv.Solve(rhs1, x);
            watch1.Stop();
            var CholTime1 = watch1.ElapsedMilliseconds;

            double Chol_nrm1 = Mat.MaxNorm(rhs1, x);

            //warm start time
            var watch2 = System.Diagnostics.Stopwatch.StartNew();
            CholSlv.Solve(rhs2, x);
            watch2.Stop();
            var CholTime2 = watch2.ElapsedMilliseconds;

            double Chol_nrm2 = Mat.MaxNorm(rhs2, x);


            // --------------- UMFPACK solver
            Solver UMFSlv = new Solver(Mat, Solver.UMFPACKSOLVER);
            //cold start time
            var watch3 = System.Diagnostics.Stopwatch.StartNew();
            UMFSlv.Solve(rhs1, x);
            watch3.Stop();
            var UMFTime1 = watch3.ElapsedMilliseconds;
            double UMF_nrm1 = Mat.MaxNorm(rhs1, x);

            //warm start time
            var watch4 = System.Diagnostics.Stopwatch.StartNew();
            UMFSlv.Solve(rhs2, x);
            watch4.Stop();
            var UMFTime2 = watch4.ElapsedMilliseconds;
            double UMF_nrm2 = Mat.MaxNorm(rhs2, x);


            // --------------- QR solver
            Solver QRSlv = new Solver(Mat, Solver.SPQRSOLVER);
            //cold start time
            var watch5 = System.Diagnostics.Stopwatch.StartNew();
            QRSlv.Solve(rhs1, x);
            watch5.Stop();
            var QRTime1 = watch5.ElapsedMilliseconds;
            double QR_nrm1 = Mat.MaxNorm(rhs1, x);

            //warm start time
            var watch6 = System.Diagnostics.Stopwatch.StartNew();
            QRSlv.Solve(rhs2, x);
            watch6.Stop();
            var QRTime2 = watch6.ElapsedMilliseconds;
            double QR_nrm2 = Mat.MaxNorm(rhs2, x);

            // --------------- Dense solver
            Solver LUSlv = new Solver(Mat, Solver.DENSESOLVER);
            //cold start time
            var watch7 = System.Diagnostics.Stopwatch.StartNew();
            LUSlv.Solve(rhs1, x);
            watch7.Stop();
            var LUTime1 = watch7.ElapsedMilliseconds;
            double LU_nrm1 = Mat.MaxNorm(rhs1, x);

            //warm start time
            var watch8 = System.Diagnostics.Stopwatch.StartNew();
            LUSlv.Solve(rhs2, x);
            watch8.Stop();
            var LUTime2 = watch8.ElapsedMilliseconds;
            double LU_nrm2 = Mat.MaxNorm(rhs2, x);


            Console.WriteLine("Test results:");
            Console.WriteLine("\tCholmod: Cold Start {0} ms, Warm Start {1} ms, Norm {2}", CholTime1, CholTime2, Chol_nrm1);
            Console.WriteLine("\tUMFPACK: Cold Start {0} ms, Warm Start {1} ms, Norm {2}", UMFTime1, UMFTime2, UMF_nrm1);
            Console.WriteLine("\tQRSOLVER: Cold Start {0} ms, Warm Start {1} ms, Norm {2}", QRTime1, QRTime2, QR_nrm1);
            Console.WriteLine("\tDENSESOLVER: Cold Start {0} ms, Warm Start {1} ms, Norm {2}", LUTime1, LUTime2, LU_nrm1);
            
            
            //return warnings
            Solver.NoWarnings = false;
        }
        
        static void Main(string[] args)
        {
            Console.WriteLine("-- Diagonal matrix");
            test1();
            
            Console.WriteLine("-- Cholezky decomposition for square unsymmetric matrix");
            test2();
            
            Console.WriteLine("-- LU decomposition for square unsymmetric matrix");
            test3();
            
            Console.WriteLine("-- QR decomposition for square unsymmteric matrix");
            test4();
            
            Console.WriteLine("-- QR for non-square matrix");
            test5();
            
            Console.WriteLine("-- Symmetric (upper triangle) matrix routines");
            test6();
            
            Console.WriteLine("-- Symmetric matrix solution with band stencil. Methods compare");
            test7(100000, new int[] { 1, 3, 5 }, new Tuple<double, double>(1, 1000));
            
            Console.WriteLine("-- Unsymmetric matrix with band stencil. Methods compare");
            test9(100000, new int[] { 1, 3, 5, 10 }, new Tuple<double, double>(1, 1000));
            
            Console.WriteLine("-- Matrix with random stencil. Methods compare");
            test8(2000, 5, new Tuple<double, double>(1,1000));
            
            Console.WriteLine("-- Random dense Matrix solution. Methods compare");
            test10(1000, new Tuple<double, double>(1, 1000));

            Console.WriteLine("Press any key ...");
            Console.ReadLine();
        }
    }
}
