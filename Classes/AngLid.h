#ifndef LIDDRIVENCAVITY_H
#define LIDDRIVENCAVITY_H
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
using namespace std;

#include "cblas.h"

#define F77NAME(x) x##_
extern "C" {
    // Performs LU factorisation of general banded matrix
    void F77NAME(dgbtrf)(const int& N, const int& M, const int& KL,
                        const int& KU, double* AB, const int& LDAB,
                        int* IPIV, int* INFO);

    // Solves pre-factored system of equations
    void F77NAME(dgbtrs)(const char& TRANS, const int& N, const int& KL,
                        const int& KU,
                        const int& NRHS, double* AB, const int& LDAB,
                        int* IPIV, double* B, const int& LDB, int* INFO);
}

class LidDrivenCavity
{
    
public:
    LidDrivenCavity();
    ~LidDrivenCavity();
    // Access and change private data 
    void SetDomainSize(double xlen, double ylen);
    void SetGridSize(int nx, int ny);
    void SetPartitionSize(int px,int py);
    void SetTimeStep(double deltat);
    void SetFinalTime(double finalt);
    void SetReynoldsNumber(double Re);
    
    //Set vectors and matrices 
    void SetBCdata(); //
    void SetConstMatrices();
    void SetVectors();
    void SetCoeffsandsizes();
    
    //Set parameters for t+dt 
    
    //Set Poisson parameters
    void PrintMatrix(int Ny, int Nx, double* H);   //Print Matrix function
    
    
    //___________________________________________________________________________________________________________________________
    //Add boundary conditions
    void TopBC();
    void BottomBC();
    void LeftBC();
    void RightBC();
    //Add boundary conditions as vectors for each case, top, bottom, left, right, on the inner vorticity.
    void AddBC();
    //___________________________________________________________________________________________________________________________
    //Transform Row to Column and opposite functions needed for the vorticity calculations at time t+dt
    void TransColtoRow(double *A, int Ny, int Nx);
    void TransRowtoCol(double *A, int Ny, int Nx);
    
    //____________________________________________________________________________________________________________________________
    /*create static matrices that need to be generated only once
     * Amatrix: upper symmetric banded matrix for the vorticity calculation at time t : omega = A * psi
     * AmatrixGB: full general banded matrix , same as Amatrix but presented fully for the LAPACK library to solve the system in the Poisson solver
     * A_Nj: used in the first order central difference, when j was varying.
     * A_Ni: used in the first order central difference, when i was varying.
     * */
    void Amatrix(); 
    void AmatrixGB();
    void A_Nj_matrix();
    void A_Ni_matrix();
    //_____________________________________________________________________________________________________________________________
    // Vorticity at time t+dt calculation: the three matrix multiplications expressed in different functions
    void Afirstmult();
    void Asecondmult();
    void Athirdmult();
    void Finalcalc();
    //______________________________________________________________________________________________________________________________
    //vorticity omega calculation at time t
    void omegat3();
    
   
    void Initialise(double xlen, double ylen,int nx, int ny,int px, int py, double deltat, double finalt, double re);
    void Integrate();


//public:
protected:
    //protected data can be accessed by the Poisson class
    double *omega_in_new = nullptr; 
    double *psi_new = nullptr; //result of the poisson
    double* A_generalbanded = nullptr;
    int lda2;
    int nsv;
    int nsv2;
    double* omega_in = nullptr; //interior points vectors
    double* psi_in = nullptr;
    int *ipiv = nullptr;
    double U;
    double dt;
    double T;
    int    Nx;
    int    Ny;
    double Lx;
    double Ly;
    double Re;
    int Px;
    int Py;
    
    //________________________________________________________________________
    double dx;  
    double dy;
    double* psi = nullptr;
    double* omega = nullptr;
   
    //omega at time t parameters
    int lda; // for the first A matrix upper symmetric
    
    double* A = nullptr; // A matrix upper symmetric  initialisation
    
    //parameters for omega at t+dt
    int lda_Nj;
    double* A_Nj = nullptr;
    double* A_Ni = nullptr;
    
    double *res_psi_j = nullptr;         //column result
    double *res_psi_i = nullptr;         //row result
    double *res_omega_j = nullptr;       //row result
    double *res_omega_i = nullptr;       //column result
    double *res_omega_in3 = nullptr;     //result of the third multiplication
     
   //multiplication results
    double *mult1 = nullptr;
    double *mult2 = nullptr;
   
    //New vectors for each BC that needs to be added to the central differrences
    double *omega_in_top = nullptr;
    double *omega_in_bottom = nullptr;
    double *omega_in_left = nullptr;
    double *omega_in_right = nullptr;


};


#endif // LIDDRIVENCAVITY_H
