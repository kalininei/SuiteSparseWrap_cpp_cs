#include <set>
#include <stdio.h>
#include <cmath>
#include <vector>
#include "cholmod_wrapper.h"
#include <ctime>

//--- verbose output
//#define INFO(S) S.Info()
//--- minimal output
#define INFO(S) 

//c procedures check
void test1(){
	int rows[] = {0, 0, 1, 1, 1, 2, 3};
	int cols[] = {1, 3, 1, 2, 3, 0, 3};
	double vals[] = {-1, 1, 4, 2, 2, 2, 1}; 
	double b[] = {1, 1, 1, 1};
	double x[] = {0, 0, 0, 0};
	void *sol;

	
	//CHOLMOD
	sol = cholmodw_init(CHOLMODSOLVER);
	cholmodw_setmat_triplets(sol, 7, 4, 4, rows, cols, vals, 0);
	cholmodw_solve(sol, b, x);
	//once more to mimic second procedure call
	cholmodw_solve(sol, b, x);
	printf("%f %f %f %f\n", x[0], x[1], x[2], x[3]);
	cholmodw_free(sol);
	
	//UMFPACK
	sol = cholmodw_init(UMFPACKSOLVER);
	cholmodw_setmat_triplets(sol, 7, 4, 4, rows, cols, vals, 0);
	cholmodw_solve(sol, b, x);
	//once more to mimic second procedure call
	cholmodw_solve(sol, b, x);
	printf("%f %f %f %f\n", x[0], x[1], x[2], x[3]);
	cholmodw_free(sol);

	//SPQR
	sol = cholmodw_init(SPQRSOLVER);
	cholmodw_setmat_triplets(sol, 7, 4, 4, rows, cols, vals, 0);
	cholmodw_solve(sol, b, x);
	//once more to mimic second procedure call
	cholmodw_solve(sol, b, x);
	printf("%f %f %f %f\n", x[0], x[1], x[2], x[3]);
	cholmodw_free(sol);

	//DENSE SOLVER
	sol = cholmodw_init(DENSESOLVER);
	cholmodw_setmat_triplets(sol, 7, 4, 4, rows, cols, vals, 0);
	cholmodw_solve(sol, b, x);
	//once more to mimic second procedure call
	cholmodw_solve(sol, b, x);
	printf("%f %f %f %f\n", x[0], x[1], x[2], x[3]);
	cholmodw_free(sol);
}

//c++ procedures check
void test2(){
	//triplets
	int rows[] = {0, 0, 1, 1, 1, 2, 3};
	int cols[] = {1, 3, 1, 2, 3, 0, 3};
	double vals[] = {-1, 1, 4, 2, 2, 2, 1}; 
	double b[] = {1, 1, 1, 1};
	double x[] = {0, 0, 0, 0};

	//rowmaj
	double all_vals[] = {0, -1, 0, 1, 0, 4, 2, 2, 2, 0, 0, 0, 0, 0, 0, 1}; 

	//csr
	int row_ptr[] = {0, 2, 5, 6, 7};
	int col_ind[] = {1, 3, 1, 2, 3, 0, 3};

	//1) CHOLMOD
	{
		SuiteSparseWrap::Solver<CHOLMODSOLVER> S1;
		S1.SetMatTriplets(7, 4, 4, rows, cols, vals, 0);
		S1.Solve(b, x);
		printf("%f %f %f %f\n", x[0], x[1], x[2], x[3]);
		S1.SetMatRowMaj(4, 4, all_vals, 0);
		S1.Solve(b, x);
		printf("%f %f %f %f\n", x[0], x[1], x[2], x[3]);
		S1.SetMatCSR(4, row_ptr, col_ind, vals, 0);
		S1.Solve(b, x);
		printf("%f %f %f %f\n", x[0], x[1], x[2], x[3]);
	}
	//2) UMFPACK
	{
		SuiteSparseWrap::Solver<UMFPACKSOLVER> S1;
		S1.SetMatTriplets(7, 4, 4, rows, cols, vals, 0);
		S1.Solve(b, x);
		printf("%f %f %f %f\n", x[0], x[1], x[2], x[3]);
		S1.SetMatRowMaj(4, 4, all_vals, 0);
		S1.Solve(b, x);
		printf("%f %f %f %f\n", x[0], x[1], x[2], x[3]);
		S1.SetMatCSR(4, row_ptr, col_ind, vals, 0);
		S1.Solve(b, x);
		printf("%f %f %f %f\n", x[0], x[1], x[2], x[3]);
	}
	//3) QRSolver
	{
		SuiteSparseWrap::Solver<SPQRSOLVER> S1;
		S1.SetMatTriplets(7, 4, 4, rows, cols, vals, 0);
		S1.Solve(b, x);
		printf("%f %f %f %f\n", x[0], x[1], x[2], x[3]);
		S1.SetMatRowMaj(4, 4, all_vals, 0);
		S1.Solve(b, x);
		printf("%f %f %f %f\n", x[0], x[1], x[2], x[3]);
		S1.SetMatCSR(4, row_ptr, col_ind, vals, 0);
		S1.Solve(b, x);
		printf("%f %f %f %f\n", x[0], x[1], x[2], x[3]);
	}
	//4) DenseSolver
	{
		SuiteSparseWrap::Solver<DENSESOLVER> S1;
		S1.SetMatTriplets(7, 4, 4, rows, cols, vals, 0);
		S1.Solve(b, x);
		printf("%f %f %f %f\n", x[0], x[1], x[2], x[3]);
		S1.SetMatRowMaj(4, 4, all_vals, 0);
		S1.Solve(b, x);
		printf("%f %f %f %f\n", x[0], x[1], x[2], x[3]);
		S1.SetMatCSR(4, row_ptr, col_ind, vals, 0);
		S1.Solve(b, x);
		printf("%f %f %f %f\n", x[0], x[1], x[2], x[3]);
	}
}

class TPseudoRandom{
	int _k;
	double fun(int k){
		double v = (double)k+1000;
		return (sin(45.11*v*v-cos((double)k)/3.0)+1)/2.0;
	}
public:
	TPseudoRandom():_k(0){}
	double GetDouble(double min, double max, int k=-1){
		if (k<0) k=_k++;
		return fun(k)*(max-min)+min;
	}
	int GetInt(int min, int max, int k=-1){
		if (k<0) k=_k++;
		return std::floor(GetDouble(min, max, k));
	}
} PseudoRandom;

double getNorm(const std::vector<int>& Ai, const std::vector<int>& Aj,
		const std::vector<double>& Ax, 
		const std::vector<double>&x, const std::vector<double>& rhs, int sym){
	if (sym!=0) { printf("getNorm for triangle matrix is not implemented.\n"); return 0; }

	std::vector<double> Res = rhs;

	for (int i=0; i<Ai.size(); ++i){
		Res[Ai[i]]-=x[Aj[i]]*Ax[i];
	}
	double ret=0;
	for (int i=0; i<Res.size(); ++i){
		if (fabs(Res[i])>ret) ret = Res[i];
	}
	return ret;
}


//big triangle good matrix
void test3(int Nr){
	printf("+++++++++++ Big triangle matrix with band stencil\n");
	int NSten = 4, Sten[] = {1,4,10,20};
	std::vector<int> Ai, Aj; std::vector<double> Ax;
	double v1, v2, v3;
	clock_t clock_start;

	//lower triangle sparse matrix build: 
	for (int i=0; i<Nr; ++i){
		for (int k=0; k<NSten; ++k){
			int j = i-Sten[k];
			if (j>=0){
				double val = PseudoRandom.GetDouble(0, 100);
				Ai.push_back(i);
				Aj.push_back(j);
				Ax.push_back(-val);
			}
		}
	}
	//diag
	std::vector<double> diag(Nr, 0);
	for (int i=0; i<Ai.size(); ++i){
		diag[Ai[i]]-=Ax[i];
		diag[Aj[i]]-=Ax[i];
	}
	for (int i=0; i<Nr; ++i){
		Ai.push_back(i);
		Aj.push_back(i);
		Ax.push_back(diag[i]+PseudoRandom.GetDouble(0,50));
	}
	std::vector<double> rhs(Nr, 1);
	std::vector<double> x(Nr, 0);

	//CHOLMOD solution
	clock_start = clock();
	SuiteSparseWrap::Solver<CHOLMODSOLVER> S;
	S.SetMatTriplets(Ax.size(), Nr, Nr, &Ai[0], &Aj[0], &Ax[0], -1);
	S.Solve(&rhs[0], &x[0]);
	double chol_time = (float(clock()-clock_start))/CLOCKS_PER_SEC;
	
	printf("CHOLMOD solution\n");
	v1=x[2];
	INFO(S);

	//QR solution
	clock_start = clock();
	SuiteSparseWrap::Solver<SPQRSOLVER> S2;
	S2.SetMatTriplets(Ax.size(), Nr, Nr, &Ai[0], &Aj[0], &Ax[0], -1);
	S2.Solve(&rhs[0], &x[0]);
	double qr_time = (float(clock()-clock_start))/CLOCKS_PER_SEC;
	
	printf("SPQR solution\n");
	v2=x[2];
	INFO(S2);
	
	//UMF solution
	clock_start = clock();
	SuiteSparseWrap::Solver<UMFPACKSOLVER> S3;
	S3.SetMatTriplets(Ax.size(), Nr, Nr, &Ai[0], &Aj[0], &Ax[0], -1);
	S3.Solve(&rhs[0], &x[0]);
	double umf_time = (float(clock()-clock_start))/CLOCKS_PER_SEC;
	
	printf("UMF solution\n");
	v3=x[2];
	INFO(S3);

	printf("Execution times:\n\tCholmod - %f sec,\n\tSPQR - %f sec,\n\tUMFPACK - %f sec.\n", 
			chol_time, qr_time, umf_time);
	if (fabs(v1-v2)>1e-6 || fabs(v2-v3)>1e-6)
		printf("TEST FAILED: different solutions for the same matrix obtained\n");
}

//big non-symmetric good matrix
void test4(int Nr){
	printf("+++++++++++ Big unsymmetric matrix with band stencil\n");
	int NSten = 4, Sten[] = {1,4,10,20};
	double v1, v2, v3;
	std::vector<int> Ai, Aj; std::vector<double> Ax;
	clock_t clock_start;

	//lower triangle sparse matrix build: 
	for (int i=0; i<Nr; ++i){
		for (int k=0; k<NSten; ++k){
			int j = i-Sten[k];
			if (j>=0){
				double val = PseudoRandom.GetDouble(0, 100);
				Ai.push_back(i); Aj.push_back(j); Ax.push_back(-val);
			}
			j = i+Sten[k];
			if (j<Nr){
				double val = PseudoRandom.GetDouble(0, 100);
				Ai.push_back(i); Aj.push_back(j); Ax.push_back(-val);
			}
		}
	}
	//diag
	std::vector<double> diag(Nr, 0);
	for (int i=0; i<Ai.size(); ++i){
		diag[Ai[i]]-=Ax[i]; diag[Aj[i]]-=Ax[i];
	}
	for (int i=0; i<Nr; ++i){
		Ai.push_back(i); Aj.push_back(i);
		Ax.push_back(diag[i]+PseudoRandom.GetDouble(0,50));
	}
	std::vector<double> rhs(Nr, 1), x(Nr, 0);

	//CHOLMOD solution
	clock_start = clock();
	SuiteSparseWrap::Solver<CHOLMODSOLVER> S;
	S.SetMatTriplets(Ax.size(), Nr, Nr, &Ai[0], &Aj[0], &Ax[0]);
	S.Solve(&rhs[0], &x[0]);
	double chol_time = (float(clock()-clock_start))/CLOCKS_PER_SEC;
	
	printf("CHOLMOD solution\n");
	v1 = x[2];
	INFO(S);

	//QR solution
	clock_start = clock();
	SuiteSparseWrap::Solver<SPQRSOLVER> S2;
	S2.SetMatTriplets(Ax.size(), Nr, Nr, &Ai[0], &Aj[0], &Ax[0]);
	S2.Solve(&rhs[0], &x[0]);
	double qr_time = (float(clock()-clock_start))/CLOCKS_PER_SEC;
	
	printf("SPQR solution\n");
	v2 = x[2];
	INFO(S2);
	
	//UMF solution
	clock_start = clock();
	SuiteSparseWrap::Solver<UMFPACKSOLVER> S3;
	S3.SetMatTriplets(Ax.size(), Nr, Nr, &Ai[0], &Aj[0], &Ax[0]);
	S3.Solve(&rhs[0], &x[0]);
	double umf_time = (float(clock()-clock_start))/CLOCKS_PER_SEC;
	
	printf("UMF solution\n");
	v3 = x[2];
	INFO(S3);

	printf("Execution times:\n\tCholmod - %f sec,\n\tSPQR - %f sec,\n\tUMFPACK - %f sec.\n", 
			chol_time, qr_time, umf_time);
	if (fabs(v1-v2)>1e-6 || fabs(v2-v3)>1e-6)
		printf("TEST FAILED: different solutions for the same matrix obtained\n");
}

//bad matrix with random stencil
void test5(int Nr){
	printf("++++++++ Matrix with random stencil\n");
	int NSten = 10;
	double v1, v2, v3, v4;
	std::vector<int> Ai, Aj; std::vector<double> Ax;
	std::vector<double> x1,x2,x3,x4;
	clock_t clock_start;

	//sparse matrix build: 
	for (int i=0; i<Nr; ++i){
		std::set<int> used;
		for (int k=0; k<NSten-1; ++k){
			int j = PseudoRandom.GetInt(0,Nr);
			if (j!=i && used.find(j) == used.end()){
				double val = PseudoRandom.GetDouble(-100, 100);
				Ai.push_back(i); Aj.push_back(j); Ax.push_back(val);
				used.insert(j);
			}
		}
		double val = PseudoRandom.GetDouble(50, 100);
		Ai.push_back(i); Aj.push_back(i); Ax.push_back(val);
	}
	std::vector<double> rhs(Nr, 1), x(Nr, 0);

	//CHOLMOD solution
	clock_start = clock();
	SuiteSparseWrap::Solver<CHOLMODSOLVER> S;
	S.SetMatTriplets(Ax.size(), Nr, Nr, &Ai[0], &Aj[0], &Ax[0]);
	S.Solve(&rhs[0], &x[0]);
	double chol_time = (float(clock()-clock_start))/CLOCKS_PER_SEC;
	
	printf("CHOLMOD solution\n");
	v1 = x[2]; x1 = x;
	INFO(S);

	//QR solution
	clock_start = clock();
	SuiteSparseWrap::Solver<SPQRSOLVER> S2;
	S2.SetMatTriplets(Ax.size(), Nr, Nr, &Ai[0], &Aj[0], &Ax[0]);
	S2.Solve(&rhs[0], &x[0]);
	double qr_time = (float(clock()-clock_start))/CLOCKS_PER_SEC;
	
	printf("SPQR solution\n");
	v2 = x[2]; x2 = x;
	INFO(S2);
	
	//UMF solution
	clock_start = clock();
	SuiteSparseWrap::Solver<UMFPACKSOLVER> S3;
	S3.SetMatTriplets(Ax.size(), Nr, Nr, &Ai[0], &Aj[0], &Ax[0]);
	S3.Solve(&rhs[0], &x[0]);
	double umf_time = (float(clock()-clock_start))/CLOCKS_PER_SEC;
	
	printf("UMF solution\n");
	v3 = x[2]; x3 = x;
	INFO(S3);

	//DENSE solution
	clock_start = clock();
	SuiteSparseWrap::Solver<DENSESOLVER> S4;
	S4.SetMatTriplets(Ax.size(), Nr, Nr, &Ai[0], &Aj[0], &Ax[0]);
	S4.Solve(&rhs[0], &x[0]);
	double lu_time = (float(clock()-clock_start))/CLOCKS_PER_SEC;
	
	printf("DENSE solution\n");
	v4 = x[2]; x4 = x;
	INFO(S4);

	printf("Execution times:\n\tCholmod - %f sec,\n\tSPQR - %f sec,\n\tUMFPACK - %f sec, \n\tDense LU - %f sec.\n", 
			chol_time, qr_time, umf_time, lu_time);
	if (fabs(v1-v2)>1e-6 || fabs(v2-v3)>1e-6 || fabs(v3-v4)>1e-6){
		printf("TEST FAILED: different solutions for the same matrix obtained\n");
		printf("x[2] = %f, %f, %f %f\n", v1, v2, v3, v4);
		printf("norms: %f, %f, %f %f\n", 
			getNorm(Ai, Aj, Ax, x1, rhs, 0),
			getNorm(Ai, Aj, Ax, x2, rhs, 0),
			getNorm(Ai, Aj, Ax, x4, rhs, 0),
			getNorm(Ai, Aj, Ax, x3, rhs, 0));
	}
}

int main(){
	test1();
	test2();
	test3(100000);
	test4(300000);
	test5(2000);

	return 0;
};







