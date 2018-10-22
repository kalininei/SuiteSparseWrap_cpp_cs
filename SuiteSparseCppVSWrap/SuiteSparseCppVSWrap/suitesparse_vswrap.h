#ifndef SUITESPARSE_VSWRAP_H
#define SUITESPARSE_VSWRAP_H

#include "windows.h"
#include <stdexcept>

#define DENSESOLVER 1
#define CHOLMODSOLVER 2
#define UMFPACKSOLVER 3
#define SPQRSOLVER 4


// =========== Initialization
//
//Initializes solver.
//methods:
//    tp = DENSESOLVER:   dense chol/lu/qr factorization solver
//    tp = CHOLMODSOLVER: cholmod sparse Cholezky factorization (A*A' algorithm is used for non-trianlge A)
//    tp = UMFPACKSOLVER: umfpack multifrontal factorization
//    tp = SPQRSOLVER:    QR factorization
typedef void* (CALLBACK* cholmodw_initTp)(int tp);
extern cholmodw_initTp cholmodw_init;

// ============ procedure to set matrix for the solver
//
// matrix from triplet.
// obj - pointer to solver
// Nnz - number of nonzero elements
// Nrow, Ncol - number of columns and rows
// rows, cols, vals - arrays of [Nnz] lenght. Each non-zero element is defined by certain entry in each of these arrays
// sym - 0 - unsymmetric matrix, -1/1 - lower/upper triangle matrix
// Example
// M = |  0 -1 0 1|
//     |  0  4 2 2|
//     |  2  0 0 0|
//     |  0  0 0 1|
// will be passed by:
// 
// void* sol = cholmodw_init();
// int rows[] = {0, 0, 1, 1, 1, 2, 3};
// int cols[] = {1, 3, 1, 2, 3, 0, 3};
// double vals[] = {-1, 1, 4, 2, 2, 2, 1}; 
// cholmodw_setmat_triplets(sol, 7, 4, 4, rows, cols, vals, 0);
typedef void (CALLBACK* cholmodw_setmat_tripletsTp)(void* obj, int Nnz, int Nrow, int Ncol, int* rows, int* cols, double* vals, int sym);
extern cholmodw_setmat_tripletsTp cholmodw_setmat_triplets;

// matrix in dense rowmajor format. All zero entries will be eliminated within the procedure
// Example
// M = |  0 -1 0 1|
//     |  0  4 2 2|
//     |  2  0 0 0|
//     |  0  0 0 1|
// will be passed by:
// 
// void* sol = cholmodw_init();
// double vals[] = {0, -1, 0, 1, 0, 4, 2, 2, 2, 0, 0, 0, 0, 0, 0, 1}; 
// cholmodw_setmat_rowmaj(sol, 4, 4, vals, 0);
typedef void (CALLBACK* cholmodw_setmat_rowmajTp)(void* obj, int Nrow, int Ncol, double* vals, int sym);
extern cholmodw_setmat_rowmajTp cholmodw_setmat_rowmaj;

//matrix in csr format
// Example
// M = |  0 -1 0 1|
//     |  0  4 2 2|
//     |  2  0 0 0|
//     |  0  0 0 1|
// will be passed by:
// 
// void* sol = cholmodw_init();
// int row_ptr[] = {0, 2, 5, 6, 7};
// int col_ind[] = {1, 3, 1, 2, 3, 0, 3};
// double vals[] = {-1, 1, 4, 2, 2, 2, 1}; 
// cholmodw_setmat_csr(sol, 4, row_ptr, col_ind, vals, 0);
typedef void (CALLBACK* cholmodw_setmat_csrTp)(void* obj, int Nrows, int* row_ptr, int* col_ind, double* vals, int sym);
extern cholmodw_setmat_csrTp cholmodw_setmat_csr;

// ==================== Solution

//solve a system Ax = rhs
//This procedure can be called multiple times with various rhs.
//Second and all further calls will be much faster, because factorization procedures will be omitted.
//Returns non zero integer if solution was obtained and 0 otherwise.
typedef int  (CALLBACK* cholmodw_solveTp)(void* obj, double* rhs, double* x);
extern cholmodw_solveTp cholmodw_solve;

// ===================== Post processing
//print solution information into stdout
typedef void (CALLBACK* cholmodw_print_infoTp)(void* obj);
extern cholmodw_print_infoTp  cholmodw_print_info;

//free all memory allocated for the solver
typedef void (CALLBACK* cholmodw_freeTp)(void* obj);
extern cholmodw_freeTp cholmodw_free;

// ===================== C++ interface
namespace SuiteSparseWrap {

	//class wrapper for procedures above
	template<int tp = SPQRSOLVER>
	class Solver {
		void* slv;
	public:
		Solver() :slv(cholmodw_init(tp)) {}
		~Solver() { cholmodw_free(slv); }

		//add matrix
		void SetMatTriplets(int Nnz, int Nrow, int Ncol, int* rows, int* cols, double* vals, int sym = 0) {
			cholmodw_setmat_triplets(slv, Nnz, Nrow, Ncol, rows, cols, vals, sym);
		}
		void SetMatRowMaj(int Nrow, int Ncol, double* vals, int sym = 0) {
			cholmodw_setmat_rowmaj(slv, Nrow, Ncol, vals, sym);
		}
		void SetMatCSR(int Nrows, int* row_ptr, int* col_ind, double* vals, int sym = 0) {
			cholmodw_setmat_csr(slv, Nrows, row_ptr, col_ind, vals, sym);
		}

		//solution
		void Solve(double* rhs, double* x) {
			if (cholmodw_solve(slv, rhs, x) == 0)
				throw std::runtime_error("Matrix solution was not obtained");
		}

		//print solution info to stdout
		void Info() {
			cholmodw_print_info(slv);
		}
	};

}//namespace

#endif
