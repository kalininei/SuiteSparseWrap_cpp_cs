#ifndef CHOLMOD_WRAPPER_H
#define CHOLMOD_WRAPPER_H

#ifdef WIN32 
	#include <stdlib.h>
	#define STDCALL __stdcall
#else
	#define STDCALL
#endif

#define DENSESOLVER 1
#define CHOLMODSOLVER 2
#define UMFPACKSOLVER 3
#define SPQRSOLVER 4

// ================= C interface
#ifdef __cplusplus
extern "C"{
#endif

// =========== Initialization
//
//Initializes solver.
//methods:
//    tp = DENSESOLVER:   dense chol/lu/qr factorization solver
//    tp = CHOLMODSOLVER: cholmod sparse Cholezky factorization (A*A' algorithm is used for non-trianlge A)
//    tp = UMFPACKSOLVER: umfpack multifrontal factorization
//    tp = SPQRSOLVER:    QR factorization
void* STDCALL cholmodw_init(int tp);

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
void STDCALL cholmodw_setmat_triplets(void* obj, int Nnz, int Nrow, int Ncol, int* rows, int* cols, double* vals, int sym);

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
void STDCALL cholmodw_setmat_rowmaj(void* obj, int Nrow, int Ncol, double* vals, int sym);

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
void STDCALL cholmodw_setmat_csr(void* obj, int Nrows, int* row_ptr, int* col_ind, double* vals, int sym);

// ==================== Solution

//solve a system Ax = rhs
//These procedure can be called multiple times with various rhs.
//Second and all further calls will be much faster, because factorization procedures will be omitted.
//Returns non zero integer if solution was obtained and 0 otherwise.
int STDCALL cholmodw_solve(void* obj, double* rhs, double* x);

// ===================== Post processing
//print solution information into stdout
void STDCALL cholmodw_print_info(void* obj);

//free all memory allocated for the solver
void STDCALL cholmodw_free(void* obj);

#ifdef __cplusplus
}//extern C
#endif



// ===================== C++ header interface
#ifdef __cplusplus
#include <stdexcept>

namespace SuiteSparseWrap{

//class wrapper for procedures above
template<int tp = SPQRSOLVER>
class Solver{
	void* slv;
public:
	Solver():slv(cholmodw_init(tp)) {}
	~Solver() { cholmodw_free(slv); }

	//add matrix
	void SetMatTriplets(int Nnz, int Nrow, int Ncol, int* rows, int* cols, double* vals, int sym=0){
		cholmodw_setmat_triplets(slv, Nnz, Nrow, Ncol, rows, cols, vals, sym);
	}
	void SetMatRowMaj(int Nrow, int Ncol, double* vals, int sym=0){
		cholmodw_setmat_rowmaj(slv, Nrow, Ncol, vals, sym);
	}
	void SetMatCSR(int Nrows, int* row_ptr, int* col_ind, double* vals, int sym=0){
		cholmodw_setmat_csr(slv, Nrows, row_ptr, col_ind, vals, sym);
	}

	//solution
	void Solve(double* rhs, double* x){
		if (cholmodw_solve(slv, rhs, x) == 0) 
			throw std::runtime_error("Matrix solution was not obtained");
	}
	
	//print solution info to stdout
	void Info(){ 
		cholmodw_print_info(slv);
	}
};

}//namespace

#endif


#endif
