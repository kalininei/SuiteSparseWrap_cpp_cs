using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Runtime.InteropServices;

namespace SuiteSparseCSWrap
{
	internal static class NativeMethods
	{
		[DllImport("cholmod_wrapper.dll")]
		internal static extern IntPtr cholmodw_init(int tp);
		[DllImport("cholmod_wrapper.dll")]
		internal static extern void cholmodw_free(IntPtr obj);
		[DllImport("cholmod_wrapper.dll")]
		internal static extern void cholmodw_setmat_triplets(
				IntPtr obj, int Nnz, int Nrow, int Ncol, int[] rows, int[] cols, double[] vals, int sym);
		[DllImport("cholmod_wrapper.dll")]
		internal static extern int cholmodw_solve(IntPtr obj, double[] rhs, double[] x);
		[DllImport("cholmod_wrapper.dll")]
		internal static extern void cholmodw_print_info(IntPtr obj);
	}

    /// <summary>
    /// Handle to a unmanaged solver objecct
    /// </summary>
    internal class CSolverHandle 
    {
        internal IntPtr Handle { get; set; }

        internal CSolverHandle(IntPtr handle) { this.Handle = handle; }

        internal bool IsInvalid
        {
            get { return (this.Handle == IntPtr.Zero); }
        }

        internal void CloseHandle()
        {
            // Do nothing if the handle is invalid
            if (this.IsInvalid) return;

            // Close the handle, logging but otherwise ignoring errors
            NativeMethods.cholmodw_free(this.Handle);

            // Set the handle to an invalid value
            this.Handle = IntPtr.Zero;
        }
    }


	/// <summary>
	/// Sparse Matrix Solver class
	/// </summary>
	public class Solver
	{
		/// <summary>
		/// set this to true to supress all warnings
		/// </summary>
		public static bool NoWarnings = false;
		
		/// <summary>
		/// Dense LU matrix solver (uses lapack procedures)
		/// </summary>
		public const int DENSESOLVER = 1;
		/// <summary>
		/// Sparse matrix cholmod solver. For symmetric positive defined matrices. 
		/// If matrix is non-symmetric implements A*A' algorithm (It's better to use UMFPACKSOLVER in such cases).
		/// </summary>
		public const int CHOLMODSOLVER = 2;
		/// <summary>
		/// Sparse matrix umfpack solver. LU factorization of non-symmetric square matrices 
		/// </summary>
		public const int UMFPACKSOLVER = 3;

		/// <summary>
		/// Sparse matrix qr factorization solver. QR factorization of arbitrary (including non-symmetric, non-square) matrices.
		/// </summary>
		public const int SPQRSOLVER = 4;

		/// <summary>
		/// Initializes the solver
		/// </summary>
		/// <param name="Mat">Left Hand Side matrix</param>
		/// <param name="tp">DENSESOLVER, CHOLMODSOLVER, UMFPACKSOLVER, SPQRSOLVER (default)</param>
		public Solver(SparseMatrix Mat, int tp = SPQRSOLVER)
		{
			c = new CSolverHandle(NativeMethods.cholmodw_init(tp));
			Tuple<int, int> dim = Mat.GetDim();
			int Nnz = Mat.GetNnz();
			
			//Show warning message if matrix is not sparse
			if (!NoWarnings && tp != DENSESOLVER && Nnz > 100 && Nnz > 0.7 * dim.Item1 * dim.Item2)
			{
				System.Windows.Forms.MessageBox.Show
					("Matrix seems to be dense. You should consider using dense matrix solver", "SuiteSparseCSWrap.dll warning");
			}

			//fill matrix
			var vals = Mat.GetAllValues();
			int[] rv = new int[Nnz];
			int[] cv = new int[Nnz];
			double[] vv = new double[Nnz];
			int ind = 0;
			foreach (var e in vals)
			{
				rv[ind] = e.Item1; cv[ind] = e.Item2; vv[ind] = e.Item3;
				++ind;
			}
			NativeMethods.cholmodw_setmat_triplets(c.Handle, Nnz, dim.Item1, dim.Item2, rv, cv, vv, Mat.GetSymType());
		}
		/// <summary>
		/// Solves a matrix for defined rhs vector
		/// </summary>
		/// <param name="rhs">Right hand side dense vector</param>
		/// <param name="x">Solution vector. Should be allocated with number of rows size</param>
		public void Solve(double[] rhs, double[] x)
		{
			int res = NativeMethods.cholmodw_solve(c.Handle, rhs, x);
			if (res == 0) throw new Exception("Error during the SLAE solution");
		}

		/// <summary>
		/// Prints solver info to stdout
		/// </summary>
		public void PrintInfo()
		{
			NativeMethods.cholmodw_print_info(c.Handle);
		}

        /// <summary>
        /// Implicit unmanaged memory clean.
        /// </summary>
        public void Dispose()
        {
            c.CloseHandle();
        }

        /// <summary>
        /// Destructor clears unmanaged memory
        /// </summary>
		~Solver()
		{
            Dispose();
		}

		//SuiteSparce Internal Data
		//cholmod_common
		private CSolverHandle c;
	}


}
