#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <fstream>
using namespace std;
#include <string>   
#define _USE_MATH_DEFINES //allows to use pi as M_PI using math,h 
#include <math.h>
#include <cblas.h>
#include <cmath>
#include "cstring"//to initalise evrrything to 0

#include "LidDrivenCavity.h"
//print matrix
#define F77NAME(x) x##_
extern "C" {
    //LU factorisation of general banded matrix
    void F77NAME(dgbtrf)(const int& N, const int& M, const int& KL, 
                        const int& KU, double* AB, const int& LDAB,
                        int* IPIV, int* INFO);

    // Solves pre-factored system of equations
    void F77NAME(dgbtrs)(const char& TRANS, const int& N, const int& KL, 
                        const int& KU,
                        const int& NRHS, double* AB, const int& LDAB,
                        int* IPIV, double* B, const int& LDB, int* INFO);
}


void FillSBMatrix(int nsv, double* A, int ldA,int Nx,int Ny, const double one_dx2,
const double one_dy2,const double two_dxdy);

void FillGBMatrix(int nsv, double* A_gb, int ldAgb,int Nx,int Ny, const double one_dx2,
const double one_dy2,const double two_dxdy);


 void LUfactorisation(double* A, int nsv2, int lda2,int KL,int KU){     
    int info;
    int *ipiv = new int[nsv2];
    
    F77NAME(dgbtrf)(nsv2, nsv2, KL, KU, A, lda2, ipiv, &info);

    if (info) {
        cout << "Failed to LU factorise matrix" << endl;
    }
}

void SolveLapack(int nsv, double* A, double* omega, double *psi, double* u, int nsv2, int lda2, int* ipiv) {
    int info;
    int nrhs = 1;
    int KL = floor(lda2/2);         // Number of lower diagonals
    int KU = floor(lda2/2);         // Number of upper diagonals

    cblas_dcopy(nsv, omega, 1, u, 1);

    F77NAME(dgbtrs)(lda2, nsv2, KL, KU, nrhs, A, 3*KL+1, ipiv, u, nsv2, &info);
    if (info) {
        cout << "Error in solve: " << info << endl;
    }
}


void print_matrix (const double *A,const int row, const int col) {
    for (int i=0;i<row;++i) {
        for (int j=0;j<col;++j){
            cout<<setw(4)<<A[i+row*j]<<"  ";   
        }
        cout<<endl;
    }
    cout<<endl;
}

int main () {
    double U=1.0;
    double Lx=1.0;
    double Ly=1.0;
    double Re=100;
    int Nx=5;//columns
    int Ny=5;//rows
    double dt=1;//time step 
    double T=100;
    double dx=Lx/(Nx-1.0);
    double dy=Ly/(Ny-1.0);

    //stream function initialisation 
    double* s=new double[(Nx)*(Ny)]; 
    //memset(s, 0, (Nx)*(Ny)*sizeof(s[0]));//to fill array with 0s - check it 
    //initialise stream function once

    for (int i=0;i<Nx*Ny;++i) {
      s[i]=0.0;
    }
    
    for (int i=1;i<Nx-1;++i) {
        for (int j=1;j<Ny-1;++j){
            s[i*Ny+j]=1.0;
        }
    }
    print_matrix(s,Ny,Nx);

    /*
    boundary consitions will be only used in Poisson solver RHS 
    keep convention of handout i.e i=columns, j=rows 
    */
    //vorticity. boundary conditions @ corners v always 0
    double* v=new double[(Nx)*(Ny)];
    /*top BC
    needs access to N_y-2 row of s-> rows*index + no. positions to last element
    */
    
    /*bottom BC
    now we need 2nd row of s -> rows*index+1 
    */
    for (int i=0;i<Nx;++i){//Nx no of columns 
        v[Ny*i+Ny-1]=-s[(Ny*i)+Ny-2]*2/dy/dy-2*U/dy;//last row of s
        v[(Ny*i)]=-s[(Ny*i)+1]*2/dy/dy;//first row of v
    }

    /**
    left needs access to 2nd (left) column 
    right needs access last to final column (Nx-2) s-> rows*columns 
    */

    for (int i=0;i<Ny;++i){
        v[i]=-s[Ny+i]*2/dx/dx;//first column of v (left)
        v[(Nx-1)*Ny+i]=-s[(Nx-2)*Ny+i]*2/dx/dx;//last column of v (right)

    }
    //produce the symmetric banded matrix A remember ldA*Number of columns 
    //ldA = Kl + Ku + 1 for symm banded - allowed in blas 
    //exclude padding rows 
    //no. of columns = [Ny-2]*[Nx-2](leading dimension of banded matrix)
    int ldA=3+(Ny-4);//0 diagonals is ny-4
    int nsv=(Ny-2)*(Nx-2);//inner stream/vorticity size 
    double *A=new double [ldA*nsv];
    const double one_dx2=-1/dx/dx;//1st KU 
    const double one_dy2=-1/dy/dy;//3rd KU
    const double two_dxdy=2/dx/dx+2/dy/dy;//diag 
    //v=A s : s is a vector of vector columns
    //configure A elements up to position 3 
    //diagonal
    FillSBMatrix(nsv,A,ldA,Nx,Ny,one_dx2,one_dy2,two_dxdy);
    print_matrix(A,ldA,nsv);

    int KL=2+(Ny-4);
    int KU=2+(Ny-4);
    int ldAgb=1+2*KL+KU;//2*Kl+Ku+1
    double *A_gb=new double [ldAgb*nsv];
    FillGBMatrix(nsv,A_gb,ldAgb,Nx,Ny,one_dx2,one_dy2,two_dxdy);
    print_matrix(A_gb,ldAgb,nsv);
    LUfactorisation(A_gb,nsv,ldAgb,KL,KU);
    
    print_matrix(A_gb,ldAgb,nsv);
     
     //print_matrix(A_gb,ldAgb,nsv);
    //solve for inner vorticity 
    //use inner columns of s in vector form 
    //operate on submatrices of s and v 
    //m*n  = Ny-2*Nx-2 (submatrix)
    double *v_inner=new double [(Ny-2)*(Nx-2)];//no need for preallocation 
    double *s_inner=new double [(Ny-2)*(Nx-2)];
   
    //initialise s_inner 
    for (int i=0;i<(Nx-1)*(Ny-1);++i) {
        s_inner[i]=1;
   
    }
    cblas_dsbmv(CblasColMajor, CblasUpper,nsv,KL,1.0,A,ldA,
                    s_inner,1,0.0,v_inner,1); // v_inner <= A s_inner + v_inner

    //Calculation of inner vorticity. After that make v by compining boundary conditions and inner
    //indices from 1 to N-1 as we only want the inner vorticities
    
    print_matrix(v_inner,Ny-2,Nx-2);
    //for loop to check if correct
    
    for (int i=1;i<Nx-1;++i) {//columns 
        for (int j=1;j<Ny-1;++j){//rows
            v[i*Ny+j]=((one_dx2)*(s[(i+1)*Ny+j]-2*s[i*Ny+j]+s[(i-1)*Ny+j])+
            ((one_dy2)*(s[i*Ny+j+1]-2*s[i*Ny+j]+s[i*Ny+j-1])));
        }
    }

    print_matrix(v,Ny,Nx);
    //Banded matrix for central difference scheme
    double *B_col=new double[3*Nx-2];
    //commplete vorticity matrix with boundaries 
    delete [] v,s,v_inner,s_inner,A,A_gb,B_col;
}

void FillGBMatrix(int nsv, double* A_gb, int ldAgb,int Nx,int Ny, const double one_dx2,
const double one_dy2,const double two_dxdy){
        
        for (int i = 0; i < nsv; ++i) {
        A_gb[i*ldAgb+(Ny-2)]=one_dx2;//doesnt change - final Kl row 
            
        for (int j=1;j<(Ny-3);++j) {//populate all diagonals that have 0s
            A_gb[i*ldAgb+j+(Ny-2)]=0.0;
        }

        if (i%(Ny-2)==0) {//every 3rd position 0 
            A_gb[i*ldAgb+(Ny-4)+1+(Ny-2)]=0.0;
        }

        else {
            A_gb[i*ldAgb+(Ny-4)+1+(Ny-2)]=one_dy2;
        }

        A_gb[i*ldAgb+(Ny-4)+2+(Ny-2)]=two_dxdy;//row of 4s
        
        if ((i+1)%(Ny-2)==0) {//every third position 0, offset 1 
            A_gb[i*ldAgb+(Ny-4)+3+(Ny-2)]=0.0;
        }
        
        else {
            A_gb[i*ldAgb+(Ny-4)+3+(Ny-2)]=one_dy2;
        }

        for (int j=1;j<(Ny-3);++j) {//populate all diagonals that have 0s
            A_gb[i*ldAgb+(Ny-4)+3+j+(Ny-2)]=0.0;
        }
        
        A_gb[i*ldAgb+(2*(Ny-4)+3)+1+(Ny-2)]=one_dx2;
    }
}


void FillSBMatrix(int nsv, double* A, int ldA,int Nx,int Ny, const double one_dx2, const double one_dy2,const  double two_dxdy){
     for (int i = 0; i < nsv; ++i) {
        A[i*ldA]=one_dx2;//doesnt change - final Kl row 
            
        for (int j=1;j<(Ny-3);++j) {//populate all diagonals that have 0s
            A[i*ldA+j]=0.0;
        }
        if (i%(Ny-2)==0) {
            A[i*ldA+(Ny-4)+1]=0.0;
         }
        else {
            A[i*ldA+(Ny-4)+1]=one_dy2;
        }
          A[i*ldA+(Ny-4)+2]=two_dxdy;//row of 4s
        }
}
