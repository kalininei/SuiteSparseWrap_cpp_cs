#include <vector>
#include <map>
#include "cholmod.h"
#include "umfpack.h"
#include "SuiteSparseQR_C.h"
#include "suitesparse_wrap.h"

//lapack procedures: lu factorization and solution
extern "C" void dgetrf_(void*, void*, void*, void*, void*, void*);
extern "C" void dgetrs_(void*, void*, void*, void*, void*, void*, void*, void*, void*);


//solver base class
struct solver_impl{
	solver_impl(): is_init(false){};
	virtual void setmat_triplets(int Nnz, int Nrow, int Ncol, int* rows, int* cols, double* vals, int sym) = 0;
	virtual void setmat_rowmaj(int Nrow, int Ncol, double* vals, int sym){
		//convert to triplets
		std::vector<int> r, c;
		std::vector<double> v;
		for (int i=0; i<Nrow; ++i){
			int jstart, jend;
			if (sym == 0){
				jstart = 0; jend = Ncol;
			} else if (sym > 0) {
				jstart = i; jend = Ncol;
			} else {
				jstart = 0; jend = i+1;
			}
			for (int j=jstart; j<jend; ++j){
				if (*vals != 0){
					r.push_back(i);
					c.push_back(j);
					v.push_back(*vals);
				}
				++vals;
			}
		}
		//set matrix
		setmat_triplets(v.size(), Nrow, Ncol, &r[0], &c[0], &v[0], sym);
	}
	virtual void setmat_csr(int Nrow, int* row_ptr, int* col_ind, double* vals, int sym){
		//convert to triplets
		std::vector<int> r, c;
		std::vector<double> v;
		int maxc = -1;
		for (int i=0; i<Nrow; ++i){
			for (int j=row_ptr[i]; j<row_ptr[i+1]; ++j){
				r.push_back(i);
				c.push_back(*col_ind++);
				if (c.back()>maxc) maxc = c.back();
				v.push_back(*vals++);
			}
		}
		//build matrix
		setmat_triplets(v.size(), Nrow, maxc+1, &r[0], &c[0], &v[0], sym);
	}

	int solve(double* rhs, double* x){
		if (!is_init) _initialize();
		return _solve(rhs, x);
	}

	virtual void info() = 0;
	virtual ~solver_impl(){};

protected:
	virtual void _initialize() = 0;
	virtual int _solve(double* rhs, double* x) = 0;
	bool is_init;
};

struct suite_sparse_solver: public solver_impl{
	suite_sparse_solver(): solver_impl(), A(0){
		cholmod_l_start(&common);
	}
	// ==================== set matrix routines
	virtual void setmat_triplets(int Nnz, int Nrow, int Ncol, int* rows, int* cols, double* vals, int sym){
		freemat();
		cholmod_triplet *T = cholmod_l_allocate_triplet(Nrow, Ncol, Nnz, sym, CHOLMOD_REAL, &common);
		T->nnz = 0;
		for (int i = 0; i<Nnz; ++i){
			if (sym>0 && rows[i]>cols[i]) continue;
			if (sym<0 && rows[i]<cols[i]) continue;
			((SuiteSparse_long*)T->i)[T->nnz] = rows[i];
			((SuiteSparse_long*)T->j)[T->nnz] = cols[i];
			((double*)T->x)[T->nnz] = vals[i];
			++T->nnz;
		}
		A = cholmod_l_triplet_to_sparse(T, 0, &common);
		cholmod_l_free_triplet(&T, &common);
	}

	virtual void info(){
		cholmod_l_print_common(0, &common);
	}
protected:
	void freemat(){
		is_init = false;
		if (A!=0) cholmod_l_free_sparse(&A, &common);
		A = 0;
	}
	virtual ~suite_sparse_solver(){
		freemat();
		cholmod_l_finish(&common);
	}
	
	cholmod_common common;
	cholmod_sparse *A;
};

struct cholmod_solver: public suite_sparse_solver{
	cholmod_solver(): suite_sparse_solver(), x(0), y(0), e(0), b(0), r(0), W(0), L(0){}
	~cholmod_solver(){
		freeinit();
	}
protected:
	void freeinit(){
		is_init = false;
		if (x!=0) cholmod_l_free_dense(&x, &common); x=0;
		if (b!=0) cholmod_l_free_dense(&b, &common); b=0;
		if (r!=0) cholmod_l_free_dense(&x, &common); x=0;
		if (L!=0) cholmod_l_free_factor(&L, &common); L=0;
		if (y!=0) cholmod_l_free_dense(&y, &common); y=0;
		if (e!=0) cholmod_l_free_dense(&e, &common); e=0;
		if (W!=0) cholmod_l_free_dense(&W, &common); W=0;
	}
	void _initialize(){
		freeinit();
		b = cholmod_l_zeros(A->nrow, 1, A->xtype, &common);
		L = cholmod_l_analyze (A, &common) ; /* analyze */
		cholmod_l_factorize(A, L, &common) ; /* factorize */
		x = y = e = W = 0;  //will be allocated within solve2
		is_init = true;
	}
	int _solve(double* rhs, double* v){
		std::copy(rhs,rhs+A->nrow,(double*)b->x);
		int ret = cholmod_l_solve2(CHOLMOD_A, L, b, NULL, &x, NULL, &y, &e, &common);
		//error return
		if (ret == 0) return 0;
		//if A is nonsymmetric: (A*A')x=b system was solved
		//=> y = A'x is the answer of the original system
		if (A->stype == 0){
			if (W==0) W = cholmod_l_allocate_dense(A->ncol, 1, A->ncol, A->xtype, &common) ;
			double zero[] = {0, 0};
			double one[]= {1, 0};
			cholmod_l_sdmult(A, 2, one, zero, x, W, &common);
			std::copy((double*)W->x, (double*)W->x+A->nrow, v);
			cholmod_l_free_dense(&W, &common);
		} else {
			std::copy((double*)x->x, (double*)x->x+A->nrow, v);
		}
		//return error if v = nan
		if (v[0]!=v[0]) return 0;
		return ret;
	}
	cholmod_dense *x, *y, *e, *b, *r, *W;
	cholmod_factor *L;
};

struct umfpack_solver: public solver_impl{
	umfpack_solver(): solver_impl(), Symbolic(0), Numeric(0) {
		umfpack_di_defaults(Control);
		Control[UMFPACK_PRL] = 2;  //verbose output
		//Control [UMFPACK_ORDERING] = UMFPACK_ORDERING_METIS;
	}
	~umfpack_solver(){
		freeinit();
	}
	struct ijv_ordered: public std::map<int, std::map<int, double> >{
		void set_val(int i, int j, double v){
			iterator fnd = this->find(j);
			if (fnd == this->end()){
				std::map<int, double> m; m[i] = v;
				this->insert(std::make_pair(j, m));
			} else {
				fnd->second[i] = v;
			}
		}
	};
	void setmat_triplets(int Nnz, int nrow, int ncol, int* rows, int* cols, double* vals, int sym){
		freeinit();
		Ap.clear(); Ai.clear(); Ax.clear();
		Nrow = nrow; Ncol = ncol;
		//fill with ordering in column major format
		ijv_ordered colorder;
		for (int i=0; i<Nnz; ++i){
			if (sym<0){
				if (rows[i]<cols[i]) continue;
				colorder.set_val(cols[i], rows[i], vals[i]);
			} else if (sym>0){
				if (rows[i]>cols[i]) continue;
				colorder.set_val(cols[i], rows[i], vals[i]);
			}
			colorder.set_val(rows[i], cols[i], vals[i]);
		}
		//fill matrix
		int k = 0;
		for (ijv_ordered::iterator it = colorder.begin(); it!=colorder.end(); ++it){
			Ap.push_back(k);
			for (std::map<int, double>::iterator it2 = it->second.begin(); it2!=it->second.end(); ++it2){
				Ai.push_back(it2->first); Ax.push_back(it2->second);
				++k;
			}
		}
		Ap.push_back(k);
	}
	void info(){
		umfpack_di_report_info(Control, Info);
	}
protected:
	void freeinit(){
		is_init = false;
		if (Symbolic!=0) umfpack_di_free_symbolic(&Symbolic); Symbolic = 0;
		if (Numeric!=0) umfpack_di_free_numeric(&Numeric); Numeric = 0;
	}
	void _initialize(){
		freeinit();
		status = umfpack_di_symbolic (Nrow, Ncol, &Ap[0], &Ai[0], &Ax[0], &Symbolic, Control, Info);
		if (status == UMFPACK_OK) 
			status = umfpack_di_numeric(&Ap[0], &Ai[0], &Ax[0], Symbolic, &Numeric, Control, Info);
		is_init = true;
	}
	int _solve(double* rhs, double* v){
		if (status == UMFPACK_OK){
			status = umfpack_di_solve(UMFPACK_A, &Ap[0], &Ai[0], &Ax[0], v, rhs, Numeric, Control, Info);
			return (status == UMFPACK_OK) ? 1 : 0;
		} else return 0;
	}
	
	void *Symbolic, *Numeric;
	std::vector<int> Ap, Ai;
	std::vector<double> Ax;
	int status;
	int Nrow, Ncol;
	double Control[UMFPACK_CONTROL];
	double Info[UMFPACK_INFO];
};


struct qr_solver: public suite_sparse_solver{
	qr_solver(): suite_sparse_solver(), x(0), b(0), y(0), QR(0){}
	~qr_solver(){
		freeinit();
	}
	//make non-symmetric matrix
	virtual void setmat_triplets(int Nnz, int Nrow, int Ncol, int* rows, int* cols, double* vals, int sym){
		freeinit();
		if (sym==0) suite_sparse_solver::setmat_triplets(Nnz, Nrow, Ncol, rows, cols, vals, sym);
		else {
			//force non-symmetric
			std::vector<int> new_rows(rows, rows+Nnz);
			std::vector<int> new_cols(cols, cols+Nnz);
			std::vector<double> new_vals(vals, vals+Nnz);
			for (int i = 0; i< Nnz; ++i){
				if (rows[i]  == cols[i]) continue;
				new_rows.push_back(cols[i]);
				new_cols.push_back(rows[i]);
				new_vals.push_back(vals[i]);
			}
			suite_sparse_solver::setmat_triplets(new_vals.size(), 
					Nrow, Ncol, &new_rows[0], &new_cols[0], &new_vals[0], 0);
		}
	}
protected:
	void freeinit(){
		is_init = false;
		if (x!=0) cholmod_l_free_dense(&x, &common); x=0;
		if (b!=0) cholmod_l_free_dense(&b, &common); b=0;
		if (y!=0) cholmod_l_free_dense(&y, &common); y=0;
		if (QR!=0) SuiteSparseQR_C_free(&QR, &common); QR=0;
	}
	void _initialize(){
		freeinit();
		b = cholmod_l_zeros(A->nrow, 1, A->xtype, &common);
		QR = SuiteSparseQR_C_factorize(SPQR_ORDERING_DEFAULT, SPQR_DEFAULT_TOL, A, &common);
		x = y = 0;  //will be allocated within solve
		is_init = true;
	}

	int _solve(double* rhs, double* v){
		if (QR == 0) return 0;
		std::copy(rhs,rhs+A->nrow,(double*)b->x);
		
		// Y = Q'*B
		y = SuiteSparseQR_C_qmult(SPQR_QTX, QR, b, &common) ;
		if (y == 0) return 0;

		// X = R\(E*Y)
		x = SuiteSparseQR_C_solve(SPQR_RETX_EQUALS_B, QR, y, &common) ;
		if (x == 0) return 0;

		//copy solution
		std::copy((double*)x->x, (double*)x->x+A->nrow, v);
		cholmod_l_free_dense(&y, &common); y = 0;
		cholmod_l_free_dense(&x, &common); x = 0;
		
		return 1;
	}
	cholmod_dense *x, *b, *y;
	SuiteSparseQR_C_factorization *QR;
};

struct dense_solver: public solver_impl{
	dense_solver(): solver_impl(), INFO(0){}
	~dense_solver(){}

	void setmat_triplets(int Nnz, int nrow, int ncol, int* rows, int* cols, double* vals, int sym){
		if (nrow!=ncol) throw std::runtime_error("Only square matrices for dense solver are allowed");
		is_init = false;
		Nrow = nrow; Ncol = ncol; Nrhs = 1;
		Matrix = std::vector<double>(Nrow*Ncol,0);
		IPIV = std::vector<int>(Nrow,0);
		int *r = rows, *c = cols; double *v = vals;
		for (int i = 0; i<Nnz; ++i){
			int &indi = *r++, &indj = *c++;
			double &v = *vals++;
			Matrix[Nrow*indi+indj] = v;
			if (sym!=0 && indi!=indj) Matrix[Nrow*indj+indi]=v;
		}
	}

	void info(){
		printf("Dense Matrix Solver: %ix%i\n", Nrow, Ncol);
		if (INFO == 0){ 
			printf("Successful\n");
		} else if (INFO<0){
			printf("Failed: %i argument had an illegal value\n", -INFO);
		} else {
			printf("Failed: U(%i,%i) is zero\n", INFO, INFO);
		}
	}
protected:
	std::vector<double> Matrix;
	std::vector<int> IPIV;
	int Nrow, Ncol, Nrhs;
	int INFO;

	void _initialize(){;
		dgetrf_(&Nrow, &Ncol, &Matrix[0], &Nrow, &IPIV[0], &INFO);
		is_init = true;
	}
	int _solve(double* rhs, double* x){
		if (INFO!=0) return 0;
		char N = 'T'; //transpose before solution since lapack has fortran matrix ordering
		std::copy(rhs, rhs+Nrow, x);
		dgetrs_(&N, &Nrow, &Nrhs, &Matrix[0], &Ncol, &IPIV[0], x, &Nrow, &INFO);
		return (INFO == 0) ? 1 : 0;
	}
};


// --------------------------- exported functions
void* STDCALL cholmodw_init(int tp){
	if (tp == CHOLMODSOLVER) return new cholmod_solver();
	else if (tp == UMFPACKSOLVER) return new umfpack_solver();
	else if (tp == SPQRSOLVER) return new qr_solver();
	else if (tp == DENSESOLVER) return new dense_solver();
	else return 0;
}

void STDCALL cholmodw_setmat_triplets(void* obj, int Nnz, int Nrow, int Ncol, int* rows, int* cols, double* vals, int sym){
	((solver_impl*)obj)->setmat_triplets(Nnz, Nrow, Ncol, rows, cols, vals, sym);
}
void STDCALL cholmodw_setmat_rowmaj(void* obj, int Nrow, int Ncol, double* vals, int sym){
	((solver_impl*)obj)->setmat_rowmaj(Nrow, Ncol, vals, sym);
}

void STDCALL cholmodw_setmat_csr(void* obj, int Nrow, int* row_ptr, int* col_ind, double* vals, int sym){
	((solver_impl*)obj)->setmat_csr(Nrow, row_ptr, col_ind, vals, sym);
}

int STDCALL cholmodw_solve(void* obj, double* rhs, double* x){
	return ((solver_impl*)obj)->solve(rhs, x);
}

void STDCALL cholmodw_print_info(void* obj){
	((solver_impl*)obj)->info();
}
void STDCALL cholmodw_free(void* obj){
	delete (solver_impl*)obj;
}

