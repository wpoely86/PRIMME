#include <iostream>
#include <cstdlib>
#include <time.h>

using namespace std;

double *matrix;
// dimension of matrix
int n = 10;

void MatrixVec(void *x, void *y, int *blockSize, struct primme_params *primme);

// IMPORTANT: will not compile in C++ otherwise
extern "C" {
void dsyev_(char *jobz,char *uplo,int *n,double *A,int *lda,double *W,double *work,int *lwork,int *info);
#include "primme.h"
}

int main(void)
{
    int ret; // for error codes
    primme_params primme;
    primme_preset_method method = DYNAMIC;

    primme_initialize(&primme);
    ret = primme_set_method(method, &primme);

    if(ret)
	cerr << "Error in setting method: code " << ret << endl;

    primme.n = n;
    primme.matrixMatvec = &MatrixVec;
    primme.numEvals = n; // number of eigs to find
    primme.printLevel = 5; // prints lots of stuff
    primme.maxBasisSize = n*n;
    primme.minRestartSize = 9;

    primme_display_params(primme);

    matrix = new double [n*n];

    srand(time(0));

    // fill random symmetric matrix
    for(int i=0;i<n;i++)
	for(int j=0;j<n;j++)
	    matrix[j+i*n] = matrix[i+j*n] = 10.0/RAND_MAX*rand();

    double *evals,*evecs,*rnorms;

    evals = new double [primme.numEvals];
    evecs = new double [primme.n*(primme.numEvals+primme.maxBlockSize)];
    rnorms = new double [primme.numEvals];

    ret = dprimme(evals, evecs, rnorms, &primme);

    if(ret)
	cerr << "Error in dprimme: code " << ret << endl;

    for(int i=0;i<primme.numEvals;i++)
	cout << "Eig: " << evals[i] << "\tResidual: " << rnorms[i] << endl;

    delete [] rnorms;
    delete [] evecs;
    delete [] evals;

    primme_Free(&primme);

//    double *eigenvalues = new double [n];
//
//    //diagonalize orignal matrix:
//    char jobz = 'V';
//    char uplo = 'U';
//
//    int lwork = 3*n - 1;
//
//    double *work = new double[lwork];
//
//    int info;
//
//    dsyev_(&jobz,&uplo,&n,matrix,&n,eigenvalues,work,&lwork,&info);
//
//    if(info)
//	cerr << "Error in lapack: code " << info << endl;
//
//    for(int i=0;i<n;i++)
//	cout << "Eigs lapack: " << eigenvalues[i] << endl;
//
//    delete [] work;
//    delete [] eigenvalues;

    delete [] matrix;

    return 0;
}

void MatrixVec(void *x, void *y, int *blockSize, struct primme_params *primme)
{
    // make double out of void pointers
    double *xd= static_cast<double *> (x);
    double *yd= static_cast<double *> (y);

    for(int i=0;i<n;i++)
    {
	yd[i] = 0;

	for(int j=0;j<n;j++)
	    yd[i] += matrix[j+i*n] * xd[j];
    }
}

