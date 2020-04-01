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
#include "cstring"//to initalise evrrything to
//print matrix
void print_matrix (const double *A,const int row, const int col) {
    for (int i=0;i<row;++i) {
        for (int j=0;j<col;++j){
        //taken as normal array now
            cout<<setw(4)<<A[i+row*j]<<"  ";   
        }
        cout<<endl;
    }
    cout<<endl;
}

void print_matrix_row (const double *A,const int row, const int col) {
    for (int i=0;i<row;++i) {
        for (int j=0;j<col;++j){
        //taken as normal array now
            cout<<setw(4)<<A[i*col+j]<<"  ";   
        }
        cout<<endl;
    }
    cout<<endl;
}

void FillGBMatrix(int nsv, double* A_gb, int ldAgb,int Nx,int Ny, const double one_dx2,
const double one_dy2,const double two_dxdy){
        
        for (int i = 0; i < nsv; ++i) {
        A_gb[i*ldAgb+2*(Ny-2)]=one_dx2;//doesnt change - final Kl row 
            
        for (int j=1;j<(Ny-3);++j) {//populate all diagonals that have 0s
            A_gb[i*ldAgb+j+2*(Ny-2)]=0.0;
        }

        if (i%(Ny-2)==0) {//every 3rd position 0 
            A_gb[i*ldAgb+(Ny-4)+1+2*(Ny-2)]=0.0;
        }

        else {
            A_gb[i*ldAgb+(Ny-4)+1+2*(Ny-2)]=one_dy2;
        }

        A_gb[i*ldAgb+(Ny-4)+2+2*(Ny-2)]=two_dxdy;//row of 4s
        
        if ((i+1)%(Ny-2)==0) {//every third position 0, offset 1 
            A_gb[i*ldAgb+(Ny-4)+3+2*(Ny-2)]=0.0;
        }
        
        else {
            A_gb[i*ldAgb+(Ny-4)+3+2*(Ny-2)]=one_dy2;
        }

        for (int j=1;j<(Ny-3);++j) {//populate all diagonals that have 0s
            A_gb[i*ldAgb+(Ny-4)+3+j+2*(Ny-2)]=0.0;
        }
        
        A_gb[i*ldAgb+(2*(Ny-4)+3)+1+2*(Ny-2)]=one_dx2;
    }
}

void FillLOCMatrix(int NB, double* A_gb, int ldAgb,int Nx,int Ny, const double one_dx2,
const double one_dy2,const double two_dxdy){
        int j=1;
        for (int i = 0; i < NB; ++i) {
        A_gb[i*ldAgb+2*(Ny-2)]=one_dx2;//doesnt change - final Kl row 
            
        for (int j=1;j<(Ny-3);++j) {//populate all diagonals that have 0s
            A_gb[i*ldAgb+j+2*(Ny-2)]=0.0;
        }

        if ((i+j*NB)%(Ny-2)==0) {//every 3rd position 0 
            A_gb[i*ldAgb+(Ny-4)+1+2*(Ny-2)]=0.0;
        }

        else {
            A_gb[i*ldAgb+(Ny-4)+1+2*(Ny-2)]=one_dy2;
        }

        A_gb[i*ldAgb+(Ny-4)+2+2*(Ny-2)]=two_dxdy;//row of 4s
        
        if (((i+1)+j*NB)%(Ny-2)==0) {//every third position 0, offset 1 
            A_gb[i*ldAgb+(Ny-4)+3+2*(Ny-2)]=0.0;
        }
        
        else {
            A_gb[i*ldAgb+(Ny-4)+3+2*(Ny-2)]=one_dy2;
        }

        for (int j=1;j<(Ny-3);++j) {//populate all diagonals that have 0s
            A_gb[i*ldAgb+(Ny-4)+3+j+2*(Ny-2)]=0.0;
        }
        
        A_gb[i*ldAgb+(2*(Ny-4)+3)+1+2*(Ny-2)]=one_dx2;
    }
}
int main () {
    double U=1.0;
    double Lx=1.0;
    double Ly=1.0;
    double Re=100;
    int Nx=8;//columns
    int Ny=8;//rows
    int NB=17;
    double dt=1;//time step 
    double T=100;
    double dx=Lx/(Nx-1.0);
    double dy=Ly/(Ny-1.0);

    int nsv=(Ny-2)*(Nx-2);
    int KL=2+(Ny-4);
    int KU=2+(Ny-4);
    const double one_dx2=-1/dx/dx;//1st KU 
    const double one_dy2=-1/dy/dy;//3rd KU
    const double two_dxdy=2/dx/dx+2/dy/dy;//diag 
    int ldAgb=1+2*KL+2*KU;//2*Kl+Ku+1
    double *A_gb=new double [ldAgb*nsv];
    double *A_loc1=new double [ldAgb*NB];
    cout<<"General banded matrix A_gb: \n";
    FillGBMatrix(nsv,A_gb,ldAgb,Nx,Ny,one_dx2,one_dy2,two_dxdy);
    print_matrix(A_gb,ldAgb,nsv);
    cout<<"Local banded matrix A_loc1: \n";
    FillLOCMatrix(NB,A_loc1,ldAgb,Nx,Ny,one_dx2,one_dy2,two_dxdy);
    print_matrix(A_loc1,ldAgb,NB);
}   