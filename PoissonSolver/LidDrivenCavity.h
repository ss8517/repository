
#pragma once

#include <string>
#include <iostream>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <map>
#include <random>
#include <time.h>
#include <cmath>
#include <math.h>

#include "cblas.h"
#include "PoissonSolver.h"
#include "LidDrivenCavity.h"
#include "mpi.h"
using namespace std;

class LidDrivenCavity
{
public:

    LidDrivenCavity(int size_,int Nx,int Ny);
    ~LidDrivenCavity();

    void SetDomainSize(double xlen, double ylen);
    void SetGridSize(int nx, int ny);
    void SetTimeStep(double deltat);
    void SetFinalTime(double finalt);
    void SetReynoldsNumber(double Re);
    void SetPartitionSize (int px, int py);
    void SetConstants ();
    void BoundaryConditions();
    void FillSBMatrix(int nsv, double* A, int ldA,int Nx,int Ny, double one_dx2,
    double one_dy2,double two_dxdy);
    void InnerVorticity();
    void FillGBMatrix();
    void FillGBmatrix_1stDerCentral_col();
    void FillGBmatrix_1stDerCentral_row();
    void print_matrix (const double *A,const int row, const int col);
    void print_matrix_row (const double *A,const int row, const int col);
    void Col2Row(double* A_col,double* B_row,int row,int col);
    void Row2Col(double* A_row,double* B_col,int row,int col);

    void Initialise(double xlen, double ylen,int nx, int ny,
    int px, int py, double deltat, double finalt, double re);

protected: 
    double* A_loc;
private:
    
    int size_;//fundamental size of all arrays
    int Nx;
    int Ny;
    int KL;
    int KU;
    int ldAgb;
    int nsv;
    double one_dx2;
    double one_dy2;
    double two_dxdy;
    
    double* v_inner=nullptr;
    double* s_inner=nullptr;
    double *mult_3=nullptr;
    
    int Px;
    int Py;
    //double *A_gb=nullptr;
    //int *ipiv=nullptr;
    double dt;
    double T;
    double Lx;
    double Ly;
    double Re;
    //****************needed for Poisson
    double dx;
    double dy;

    double* v_BC_top=nullptr;
    double* v_BC_bottom=nullptr;
    double* v_BC_left=nullptr;
    double* v_BC_right=nullptr;
  
    int ldA;
    double *A=nullptr;
    
    int ldB;
    double one_2dy;
    double one_2dx;
    double *B_col=nullptr;
    double *B_row=nullptr;

    double *v_inner_col=nullptr;
    double *s_inner_col=nullptr;
    
    double *s_inner_row=nullptr;
    double *s_inner_row_result1=nullptr;
    double *s_inner_row_result2=nullptr;

    double *v_inner_row=nullptr;
    double *v_inner_row_result1=nullptr;
    double *v_inner_row_result2=nullptr;

    double *mult_1=nullptr;
    double *mult_2=nullptr;  

    double *C_sb=nullptr;
    double one_Redx2;
    double one_Redy2;
    double two_Redxdy;
    
    //variables to be given as Arguments to Poisson
    
    double* B_loc;
    double*B_loc_red=nullptr;

    int* desca;
    int* descb;

    int*    ipiv;
    int LWORK;
    double* WORK;
    double* AF;
    PoissonSolver* psolver;
    
};


