#pragma once

#include "mpi.h"
#include "PoissonSolver.h"
class PoissonSolver  {
    public:
        PoissonSolver();
        ~PoissonSolver();
        void print_matrix (const double *A,const int row, const int col);
        
        void BuildLocalMatrices (double* s_inner,double* A_loc, double* B_loc,double* B_loc_red, double* mult_3,int size_,int ldAgb,int Nx,int Ny,double one_dx2,double one_dy2,double two_dxdy,int* ipiv,int KL,int KU,
        int* desca,double* AF,double* WORK,int* descb)  ;

        void LUfactorisation(double* A_loc, int N, int* ipiv,int BWU,int BWL
        ,int* desca,double* AF, int LAF,double* WORK,int LWORK,int mycol);

        void LinearSolver (double* A_loc,double* B_loc, int N, int* ipiv,int BWU,int BWL
        ,int* desca,int* descb,double* AF, int LAF,double* WORK,int LWORK,int mycol);
};




