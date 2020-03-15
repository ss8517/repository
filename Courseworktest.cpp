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
void LUfactorisation(double* A, int nsv2, int lda2,int KL,int KU);
void SolveLapack(int nsv, double* A, double* omega, double *psi, double* u, int nsv2, int lda2, int* ipiv);
void Col2Row(double* A_col,double* B_row,int row,int col);
void Row2Col(double* A_row,double* B_col,int row,int col);

int main () {
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
    double* v=new double[(Nx)*(Ny)];//just for for loop check
    double* s=new double[(Nx)*(Ny)]; 
    memset(s, 0, (Nx)*(Ny)*sizeof(s[0]));//to fill array with 0s - check it 
    //initialise stream function once

    for (int i=1;i<Nx-1;++i) {
        for (int j=1;j<Ny-1;++j){
            s[i*Ny+j]=1.0;
        }
    }
    print_matrix(s,Ny,Nx);
    
    double *s_inner=new double [(Ny-2)*(Nx-2)];
    //initialise s_inner
    for (int i=0;i<(Nx-1)*(Ny-1);++i) {
        s_inner[i]=1;
    }

    cout<<"s_inner is: "<<endl;
    print_matrix(s_inner,Ny-2,Nx-2);

     /**
     * @brief Computes Boundary conditions in distinct vectors
     * of size (Ny-2)*(Nx-2) from s 
     * @param v_BC_top 
     * @param v_BC_bottom
     * @param v_BC_left
     * @param v_BC_right
     */

    double* v_BC_top=new double[(Nx-2)*(Ny-2)];
    double* v_BC_bottom=new double[(Nx-2)*(Ny-2)];
    double* v_BC_left=new double[(Nx-2)*(Ny-2)];
    double* v_BC_right=new double[(Nx-2)*(Ny-2)];
    
    BoundaryConditions(s_inner,v_BC_top,v_BC_bottom,v_BC_left,v_BC_right,Nx,Ny,dx,dy);
    cout<<"The Bopundary conditions are: \n";
    cout<<"Top: \n";
    print_matrix(v_BC_top,Ny-2,Nx-2);
    cout<<"Bottom: \n";
    print_matrix(v_BC_bottom,Ny-2,Nx-2);
    cout<<"Left: \n";
    print_matrix(v_BC_left,Ny-2,Nx-2);
    cout<<"Right: \n";
    print_matrix(v_BC_right,Ny-2,Nx-2);
    
    //produce the symmetric banded matrix A remember ldA*Number of columns 
    //ldA = Kl + Ku + 1 for symm banded - allowed in blas 
    //exclude padding rows 
    //no. of columns = [Ny-2]*[Nx-2](leading dimension of banded matrix)
    int KL=2+(Ny-4);
    int KU=2+(Ny-4);
    int ldA=3+(Ny-4);//0 diagonals is ny-4
    int nsv=(Ny-2)*(Nx-2);//inner stream/vorticity size 
    double *A=new double [ldA*nsv];
    const double one_dx2=-1/dx/dx;//1st KU 
    const double one_dy2=-1/dy/dy;//3rd KU
    const double two_dxdy=2/dx/dx+2/dy/dy;//diag 
    //v=A s : s is a vector of vector columns
    //configure A elements up to position 3 
    //diagonal
    double *v_inner=new double [(Ny-2)*(Nx-2)];//no need for preallocation 
    //empty arrays to store blas multiplication results
    
    FillSBMatrix(nsv,A,ldA,Nx,Ny,one_dx2,one_dy2,two_dxdy);
    InnerVorticity(nsv,KL,A,ldA,s_inner,v_inner);// v_inner <= A s_inner + v_inner
    cout<<"Inner vorticity @ t \n";
    print_matrix(v_inner,Ny-2,Nx-2);

     //for loops for check
    for (int i=1;i<Nx-1;++i) {
        for (int j=1;j<Ny-1;++j){//rows
            v[i*Ny+j]=((one_dx2)*(s[(i+1)*Ny+j]-2*s[i*Ny+j]+s[(i-1)*Ny+j])+
            ((one_dy2)*(s[i*Ny+j+1]-2*s[i*Ny+j]+s[i*Ny+j-1])));
        }
    }
    cout<<"Inner vorticity @ t (for loops):  \n";
    print_matrix(v,Ny,Nx);
    
    /**
     * @brief makes General banded matrices vorticity at t + Dt  
     * @param ldB - 1+kl+ku , number of rows, padding roes excluded
     * @param nsvB number of columns 
     * 
    */

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
    cout<<"1st mult s_inner_col: \n";
    cblas_dgbmv(CblasColMajor,CblasNoTrans,nsvB,nsvB,1,1,1.0,B_col,ldB,s_inner,1,0.0,s_inner_col,1);
    print_matrix(s_inner_col,Ny-2,Nx-2);
    
    //transpose v
    cout<<"1st mult v_inner_row_result2: \n";
    Col2Row(v_inner,v_inner_row,Nx-2,Ny-2);
    cblas_dgbmv(CblasColMajor,CblasNoTrans,nsvB,nsvB,1,1,1.0,B_row,ldB,v_inner_row,1,0.0,v_inner_row_result1,1);
    Row2Col(v_inner_row_result1,v_inner_row_result2,Ny-2,Nx-2);//v_inner_row_result2 is in column form
    print_matrix(v_inner_row_result2,Ny-2,Nx-2);
    //Add the boundary conditions

    cout<<"1st mult V_inner_row_result2 after BCs added: \n";
    cblas_daxpy((Nx-2)*(Ny-2),1.0,v_BC_left,1,v_inner_row_result2,1);//v_inner updated
    cblas_daxpy((Nx-2)*(Ny-2),-1.0,v_BC_right,1,v_inner_row_result2,1);//
    print_matrix(v_inner_row_result2,Ny-2,Nx-2);
    
    cout<<"Mult_1: s_inner_col * v_inner_row_result2 : \n";
    cblas_dsbmv(CblasColMajor, CblasUpper,nsv, 0, 1.0, s_inner_col, 1,v_inner_row_result2,1, 0.0,mult_1, 1);
    print_matrix(mult_1,Ny-2,Nx-2);

    cout<<"2nd mult s_inner_row_result2: \n";
    Col2Row(s_inner,s_inner_row,Nx-2,Ny-2);
    cblas_dgbmv(CblasColMajor,CblasNoTrans,nsvB,nsvB,1,1,1.0,B_row,ldB,s_inner_row,1,0.0,s_inner_row_result1,1);
    Row2Col(s_inner_row_result1,s_inner_row_result2,Ny-2,Nx-2);
    print_matrix(s_inner_row_result2,Ny-2,Nx-2);
    
    cout<<"2nd mult v_inner_col: \n";
    cblas_dgbmv(CblasColMajor,CblasNoTrans,nsvB,nsvB,1,1,1.0,B_col,ldB,v_inner,1,0.0,v_inner_col,1);
    print_matrix(v_inner_col,Ny-2,Nx-2);

    cout<<"2nd mult v_inner_col  after BCs added: \n";
    cblas_daxpy((Nx-2)*(Ny-2),1.0,v_BC_bottom,1,v_inner_col,1);//v_inner updated
    cblas_daxpy((Nx-2)*(Ny-2),-1.0,v_BC_top,1,v_inner_col,1);//
    print_matrix(v_inner_col,Ny-2,Nx-2);

    cout<<"Mult_2: s_inner_row_result2 * v_inner_col: \n";
    cblas_dsbmv(CblasColMajor, CblasUpper,nsv, 0, 1.0, s_inner_row_result2, 1,v_inner_col,1, 0.0,mult_2, 1);
    print_matrix(mult_2,Ny-2,Nx-2);

    /**
     * @brief Poisson solver 
     */
    
    int ldAgb=1+2*KL+KU;//2*Kl+Ku+1
    double *A_gb=new double [ldAgb*nsv];
    FillGBMatrix(nsv,A_gb,ldAgb,Nx,Ny,one_dx2,one_dy2,two_dxdy);
    LUfactorisation(A_gb,nsv,ldAgb,KL,KU);

    delete [] v_BC_top,v_BC_bottom,v_BC_left,v_BC_right;
    delete [] v,s,v_inner,s_inner,A,A_gb,B_col,B_row;
    delete [] s_inner_col,v_inner_col;
    delete [] v_inner_row,v_inner_row_result1,v_inner_row_result2;
    delete [] s_inner_row,s_inner_row_result1,s_inner_row_result2;
    delete [] mult_1,mult_2;
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
    //v[(Nx-1)*Ny+i]=-s[(Nx-2)*Ny+i]*2/dx/dx;
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

 void LUfactorisation(double* A, int nsv, int lda,int KL,int KU){     
    int info;
    int *ipiv = new int[nsv];
    
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
           B[i*ldB]   =-one_2dy;
        }
        B[i*ldB+1] = 0.0; 
        if ((i+1)%(Ny-2)==0) {//every 3rd position 0 
            B[i*ldB+2]   =0.0;
        }
        else {
           B[i*ldB+2] = one_2dy;
        }
    }
}

void FillGBmatrix_1stDerCentral_row(double *B, int ldB, int nsv,int Nx,double one_2dx){  
    for (int i=0; i<nsv; i++){
        if (i%(Nx-2)==0) {//every 3rd position 0 
            B[i*ldB]   =0.0;
        }
        else {
           B[i*ldB]   =-one_2dx;
        }
        B[i*ldB+1] = 0.0; 
        if ((i+1)%(Nx-2)==0) {//every 3rd position 0 
            B[i*ldB+2]   =0.0;
        }
        else {
           B[i*ldB+2] = one_2dx;
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

