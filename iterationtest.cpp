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

//#include "LidDrivenCavity.h"
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

void BoundaryConditions(double* s_inner, double* v_BC_top,double*v_BC_bottom,
double* v_BC_left,double* v_BC_right,int Nx,int Ny,double dx,double dy);

void FillSBMatrix(int nsv, double* A, int ldA,int Nx,int Ny, const double one_dx2,
const double one_dy2,const double two_dxdy);
void InnerVorticity(const int nsv,const int KL,double* A,const int ldA,double* s_inner,double* v_inner);
void FillGBMatrix(int nsv, double* A_gb, int ldAgb,int Nx,int Ny, const double one_dx2,
const double one_dy2,const double two_dxdy);
void FillGBmatrix_1stDerCentral_col(double *B, int ldB, int nsv,int Ny,double one_2dy);
void FillGBmatrix_1stDerCentral_row(double *B, int ldB, int nsv,int Nx,double one_2dx);
void print_matrix (const double *A,const int row, const int col);
void print_matrix_row (const double *A,const int row, const int col);
void LUfactorisation(double* A, int nsv2, int lda2,int KL,int KU,int *ipiv);
void PoissonSolver(int ldA, int nsv, double* A, double* v , double* u,int KL,int KU,int* ipiv);
void SolveLapack(int nsv, double* A, double* omega, double *psi, double* u, int nsv2, int lda2, int* ipiv);
void Col2Row(double* A_col,double* B_row,int row,int col);
void Row2Col(double* A_row,double* B_col,int row,int col);

int main () {

    double Lx=1.0;
    double Ly=1.0;
    double Re=100;
    int Nx=5;//columns
    int Ny=5;//rows
    double dt=.0001;//time step 
    double t=dt;
    double T= 30.;
    double dx=Lx/(Nx-1.0);
    double dy=Ly/(Ny-1.0);
    double *s_inner=new double [(Ny-2)*(Nx-2)];
    memset(s_inner, 0.0, (Nx-2)*(Ny-2)*sizeof(s_inner[0]));
    /*
    for (int i=0;i<(Nx-2)*(Ny-2);++i) {
        s_inner[i]=1.0;
    }
    */
    double *v_inner=new double [(Ny-2)*(Nx-2)];
    double* v_BC_top=new double[(Nx-2)*(Ny-2)];
    double* v_BC_bottom=new double[(Nx-2)*(Ny-2)];
    double* v_BC_left=new double[(Nx-2)*(Ny-2)];
    double* v_BC_right=new double[(Nx-2)*(Ny-2)];

    int KL=2+(Ny-4);
    int KU=2+(Ny-4);
    int ldA=3+(Ny-4);//0 diagonals is ny-4
    int nsv=(Ny-2)*(Nx-2);//inner stream/vorticity size 
    double *A=new double [ldA*nsv];
    const double one_dx2=-1/dx/dx;//1st KU 
    const double one_dy2=-1/dy/dy;//3rd KU
    const double two_dxdy=2/dx/dx+2/dy/dy;//diag 
    FillSBMatrix(nsv,A,ldA,Nx,Ny,one_dx2,one_dy2,two_dxdy);

    int ldB=3;
    int nsvB=(Ny-2)*(Nx-2);
    double one_2dy=1.0/2.0/dy;
    double one_2dx=1.0/2.0/dx;
    double *B_col=new double[ldB*nsvB];
    double *B_row=new double[ldB*nsvB];
    FillGBmatrix_1stDerCentral_col(B_col,ldB,nsvB,Ny,one_2dy);
    FillGBmatrix_1stDerCentral_row(B_row,ldB,nsvB,Nx,one_2dx);

    double *v_inner_col=new double [(Ny-2)*(Nx-2)];
    double *s_inner_col=new double [(Ny-2)*(Nx-2)];

    double *s_inner_row=new double [(Ny-2)*(Nx-2)];
    double *s_inner_row_result1=new double [(Ny-2)*(Nx-2)];
    double *s_inner_row_result2=new double [(Ny-2)*(Nx-2)];

    double *v_inner_row=new double [(Ny-2)*(Nx-2)];
    double *v_inner_row_result1=new double [(Ny-2)*(Nx-2)];
    double *v_inner_row_result2=new double [(Ny-2)*(Nx-2)];

    double *mult_1=new double [(Ny-2)*(Nx-2)];
    double *mult_2=new double [(Ny-2)*(Nx-2)];
    double *mult_3=new double [(Ny-2)*(Nx-2)];

    double *C_sb=new double [ldA*nsv];
    const double one_Redx2=1/Re/dx/dx;//1st KU 
    const double one_Redy2=1/Re/dy/dy;//3rd KU
    const double two_Redxdy=-2/Re/dx/dx-2/Re/dy/dy;//diag 
    FillSBMatrix(nsv, C_sb,ldA,Nx,Ny,one_Redx2,one_Redy2,two_Redxdy);

    int ldAgb=1+2*KL+KU;//2*Kl+Ku+1
    double *A_gb=new double [ldAgb*nsv];
    int *ipiv = new int[nsv];
    FillGBMatrix(nsv,A_gb,ldAgb,Nx,Ny,one_dx2,one_dy2,two_dxdy);
    LUfactorisation(A_gb,nsv,ldAgb,KL,KU,ipiv);//happens only once 

    while (t<T) {

        memset(s_inner_col, 0.0, (Nx-2)*(Ny-2)*sizeof(s_inner_col[0]));
        memset(s_inner_row_result1, 0.0, (Nx-2)*(Ny-2)*sizeof(s_inner_row_result1[0]));
        memset(s_inner_row_result2, 0.0, (Nx-2)*(Ny-2)*sizeof(s_inner_row_result2[0]));
        memset(v_inner_col, 0.0, (Nx-2)*(Ny-2)*sizeof(v_inner_col[0]));
        memset(v_inner_row_result1, 0.0, (Nx-2)*(Ny-2)*sizeof(v_inner_row_result1[0]));
        memset(v_inner_row_result2, 0.0, (Nx-2)*(Ny-2)*sizeof(v_inner_row_result2[0]));
       
        memset(mult_1, 0.0, (Nx-2)*(Ny-2)*sizeof(mult_1[0]));
        memset(mult_2, 0.0, (Nx-2)*(Ny-2)*sizeof(mult_2[0]));
        memset(mult_3, 0.0, (Nx-2)*(Ny-2)*sizeof(mult_2[0]));

        //update boundary conditions
        BoundaryConditions(s_inner,v_BC_top,v_BC_bottom,v_BC_left,v_BC_right,Nx,Ny,dx,dy);
        //Calculate inner vorticity 
        InnerVorticity(nsv,KL,A,ldA,s_inner,v_inner);// v_inner <= A s_inner + v_inner
        
        //s_inner col  major 
        cblas_dgbmv(CblasColMajor,CblasNoTrans,nsvB,nsvB,1,1,1.0,B_col,ldB,s_inner,1,0.0,s_inner_col,1);
        
        //v_inner row major 
        Col2Row(v_inner,v_inner_row,Nx-2,Ny-2);
        cblas_dgbmv(CblasColMajor,CblasNoTrans,nsvB,nsvB,1,1,1.0,B_row,ldB,v_inner_row,1,0.0,v_inner_row_result1,1);
        Row2Col(v_inner_row_result1,v_inner_row_result2,Ny-2,Nx-2);//v_inner_row_result2 is in column form
        //add BC for v_inner row major 
        cblas_daxpy((Nx-2)*(Ny-2),-one_2dx,v_BC_left,1,v_inner_row_result2,1);
        cblas_daxpy((Nx-2)*(Ny-2), one_2dx,v_BC_right,1,v_inner_row_result2,1);

        //1st mult 
        cblas_dsbmv(CblasColMajor, CblasUpper,nsv, 0, 1.0, s_inner_col, 1,v_inner_row_result2,1, 0.0,mult_1, 1);
        
        //s_inner row major 
        Col2Row(s_inner,s_inner_row,Nx-2,Ny-2);
        cblas_dgbmv(CblasColMajor,CblasNoTrans,nsvB,nsvB,1,1,1.0,B_row,ldB,s_inner_row,1,0.0,s_inner_row_result1,1);
        Row2Col(s_inner_row_result1,s_inner_row_result2,Ny-2,Nx-2);
        
        //v_inner col major 
        cblas_dgbmv(CblasColMajor,CblasNoTrans,nsvB,nsvB,1,1,1.0,B_col,ldB,v_inner,1,0.0,v_inner_col,1);
        //add BCs v_inner col major 
        cblas_daxpy((Nx-2)*(Ny-2),-one_2dy,v_BC_bottom,1,v_inner_col,1);
        cblas_daxpy((Nx-2)*(Ny-2),one_2dy,v_BC_top,1,v_inner_col,1);

        //mult 2 --- check whether mult_i have to be initalised to 0 each time at the beginning 
        cblas_dsbmv(CblasColMajor, CblasUpper,nsv, 0, 1.0, s_inner_row_result2, 1,v_inner_col,1, 0.0,mult_2, 1);

        //mult 3
        cblas_dsbmv(CblasColMajor, CblasUpper,nsv,KL,1.0,C_sb,ldA,
                    v_inner,1,0.0,mult_3,1);

        //add BC for second derivative 
        cblas_daxpy((Nx-2)*(Ny-2),one_Redy2,v_BC_bottom,1,mult_3,1);
        cblas_daxpy((Nx-2)*(Ny-2),one_Redy2,v_BC_top,1,mult_3,1);
        cblas_daxpy((Nx-2)*(Ny-2),one_Redx2,v_BC_left,1,mult_3,1);
        cblas_daxpy((Nx-2)*(Ny-2),one_Redx2,v_BC_right,1,mult_3,1);
        
        //Addition of the three mult - ouput matrix mult_3
        cblas_daxpy((Ny-2)*(Nx-2),-1.0,mult_1,1,mult_3,1);
        cblas_daxpy((Ny-2)*(Nx-2),1.0,mult_2,1,mult_3,1);

        //evaluate v @ t + dt 
        cblas_dscal((Ny-2)*(Nx-2),dt,mult_3,1);//can be placed in daxpy 
        cblas_daxpy((Ny-2)*(Nx-2),1.0,v_inner,1,mult_3,1);
        
        //solve system 
        PoissonSolver(ldAgb,nsv,A_gb,mult_3,s_inner,KL,KU,ipiv);
        //update time step 
        t+=dt;
        
    }
    
    print_matrix(s_inner,Ny-2,Nx-2);

    delete [] v_BC_top,v_BC_bottom,v_BC_left,v_BC_right;
    delete [] A,A_gb,B_col,B_row,C_sb;
    delete [] v_inner,s_inner;
    delete [] s_inner_col,v_inner_col;
    delete [] v_inner_row,v_inner_row_result1,v_inner_row_result2;
    delete [] s_inner_row,s_inner_row_result1,s_inner_row_result2;
    delete [] mult_1,mult_2,mult_3;
    delete [] ipiv;

    return 0;
}


void BoundaryConditions(double* s_inner, double* v_BC_top,double*v_BC_bottom,
double* v_BC_left,double* v_BC_right,int Nx,int Ny,double dx,double dy){
    
    double U=1.0;
    
    memset(v_BC_top, 0.0, (Nx-2)*(Ny-2)*sizeof(v_BC_top[0]));
    memset(v_BC_bottom, 0.0, (Nx-2)*(Ny-2)*sizeof(v_BC_bottom[0]));
    memset(v_BC_left, 0.0, (Nx-2)*(Ny-2)*sizeof(v_BC_left[0]));
    memset(v_BC_right, 0.0, (Nx-2)*(Ny-2)*sizeof(v_BC_right[0]));
    
    for (int i=0;i<Nx-2;++i){
        v_BC_top[(Ny-2)*i+Ny-2-1]=-s_inner[((Ny-2)*i)+Ny-3]*2/dy/dy-2*U/dy;//top
    }

    for (int i=0;i<Nx-2;++i){
        v_BC_bottom[((Ny-2)*i)]=-s_inner[(Ny-2)*i]*2/dy/dy;
    }  

    for (int i=0;i<Ny-2;++i){
        v_BC_left[i]=-s_inner[i]*2/dx/dx;
    }  


    for (int i=0;i<Ny-2;++i){
        v_BC_right[(Nx-2-1)*(Ny-2)+i]=-s_inner[(Nx-2-1)*(Ny-2)+i]*2/dx/dx;
    }
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

 void LUfactorisation(double* A, int nsv, int lda,int KL,int KU,int *ipiv){     
    int info;
   
    
    F77NAME(dgbtrf)(nsv, nsv, KL, KU, A, lda, ipiv, &info);

    if (info) {
        cout << "Failed to LU factorise matrix" << endl;
    }
}

void FillSBMatrix(int nsv, double* A, int ldA,int Nx,int Ny, const double one_dx2,
 const double one_dy2,const  double two_dxdy){
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

void FillGBmatrix_1stDerCentral_col(double *B, int ldB, int nsv,int Ny,double one_2dy){  
    for (int i=0; i<nsv; i++){
        if (i%(Ny-2)==0) {//every 3rd position 0 
            B[i*ldB]   =0.0;
        }
        else {
           B[i*ldB]   =one_2dy;
        }
        B[i*ldB+1] = 0.0; 
        if ((i+1)%(Ny-2)==0) {//every 3rd position 0 
            B[i*ldB+2]   =0.0;
        }
        else {
           B[i*ldB+2] = -one_2dy;
        }
    }
}

void FillGBmatrix_1stDerCentral_row(double *B, int ldB, int nsv,int Nx,double one_2dx){  
    for (int i=0; i<nsv; i++){
        if (i%(Nx-2)==0) {//every 3rd position 0 
            B[i*ldB]   =0.0;
        }
        else {
           B[i*ldB]   =one_2dx;
        }
        B[i*ldB+1] = 0.0; 
        if ((i+1)%(Nx-2)==0) {//every 3rd position 0 
            B[i*ldB+2]   =0.0;
        }
        else {
           B[i*ldB+2] = -one_2dx;
        }
    }
}

void print_matrix (const double *A,const int row, const int col) {
    for (int i=0;i<row;++i) {
        for (int j=0;j<col;++j){
            cout<<setw(7)<<A[i+row*j]<<"  ";   
        }
        cout<<endl;
    }
    cout<<endl;
}

void print_matrix_row (const double *A,const int row, const int col) {//print matrix that is in row major format
    for (int i=0;i<row;++i) {
        for (int j=0;j<col;++j){
        //taken as normal array now
            cout<<setw(7)<<A[i*col+j]<<"  ";   
        }
        cout<<endl;
    }
    cout<<endl;
}
void PoissonSolver(int ldA, int nsv, double* A, double* v , double* u,int KL,int KU,int* ipiv) {
    int info;
    int nrhs = 1;
    int ldu =nsv;

    cblas_dcopy(nsv, v, 1, u, 1);

    F77NAME(dgbtrs)('N', nsv, KL, KU, nrhs, A, ldA, ipiv, u, ldu, &info);
    if (info) {
        cout << "Error in solve: " << info << endl;
    }
}

void InnerVorticity(const int nsv,const int KL,double* A,const int ldA,double* s_inner,double* v_inner) {
     cblas_dsbmv(CblasColMajor, CblasUpper,nsv,KL,1.0,A,ldA,
                    s_inner,1,0.0,v_inner,1);
}

void Col2Row(double* A_col,double* B_row,int row,int col) {
    for (int i=0;i<col;++i) {
        for(int j=0;j<row;++j){
            B_row[col*j+i]=A_col[row*i+j];
        }
    }
}

void Row2Col(double* A_row,double* B_col,int row,int col) {
    for (int i=0;i<col;++i) {
        for(int j=0;j<row;++j){
            B_col[row*i+j]=A_row[col*j+i];
        }
    }
}

