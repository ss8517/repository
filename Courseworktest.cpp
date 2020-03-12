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

int main () {
    double U=1.0;
    double Lx=1.0;
    double Ly=1.0;
    double Re=100;
    int Nx=6;//columns
    int Ny=8;//rows
    double dt=1;//time step 
    double T=100;
    double dx=Lx/(Nx-1);
    double dy=Ly/(Ny-1);

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
        v[Ny*i+Ny-1]=-s[(Ny*i)+Ny-2]*2/dy/dy-2*U/dy;//last row of 
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

    //print v to check boundary conditions 
    for (int i=0;i<Ny;++i) {
        for (int j=0;j<Nx;++j){
        //taken as normal array now
            cout<<setw(4)<<v[i+Ny*j]<<"  ";   
        }
        cout<<endl;
    }

    //produce the symmetric banded matrix A remember ldA*Number of columns 
    //ldA = Kl + Ku + 1 for symm banded - allowed in blas 
    //exclude padding rows 
    //no. of columns = [Ny-2]*[Nx-2](leading dimension of banded matrix)
    int ldA=4;
    int nsv=(Ny-2)*(Nx-2);//inner stream/vorticity size 
    double *A=new double [ldA*nsv];
    const double one_dx2=-1/dx/dx;//1st KU 
    const double one_dy2=-1/dy/dy;//3rd KU
    const double two_dxdy=2/dx/dx+2/dy/dy;//diag 
    //v=A s : s is a vector of vector columns
    //configure A elements up to position 3 
    //diagonal
    A[0*ldA+3]=two_dxdy;
    A[1*ldA+3]=two_dxdy;
    A[2*ldA+3]=two_dxdy;
    //1st KU
     A[1*ldA+2]=one_dy2;
     A[2*ldA+2]=one_dy2;
     //2nd Ku
     A[2*ldA+2]=0.0;
      for (int i = 3; i < nsv; ++i) {
          A[i*ldA]=one_dx2;
          A[i*ldA+1]=0.0;
          if (i%3==0) {
              A[i*ldA+2]=0.0;
          }
          else {
               A[i*ldA+2]=one_dy2;
          }
          A[i*ldA+3]=two_dxdy;
        }
    
    //solve for inner vorticity 
    //use inner columns of s in vector form 
    //operate on submatrices of s and v 
    //m*n  = Ny-2*Nx-2 (submatrix)
    double *v_inner=new double [(Ny-2)*(Nx-2)];//no need for preallocation 
    double *s_inner=new double [(Ny-2)*(Nx-2)];
    //copy inner columns of s into s_inner 
     for (int i=1;i<Nx-1;++i) {
        for (int j=1;j<Ny-1;++j){
           s_inner[(i-1)*Ny+j-1]=s[i*Ny+j];
        }
     }


    //print s_inner to check 
      for (int i=0;i<Ny-2;++i) {
        for (int j=0;j<Nx-2;++j){
        //taken as normal array now
            cout<<setw(4)<<s_inner[i+Ny*j]<<"  ";
            
        }
        cout<<endl;
    }
    cout<<endl;
    
    cout<<endl;
    const int KL=3; 
    cblas_dsbmv(CblasColMajor, CblasUpper,nsv,KL,1.0,A,ldA,
                    s_inner,1,0.0,v_inner,1); // v_inner <= A s_inner + v_inner

    //Calculation of inner vorticity. After that make v by compining boundary conditions and inner
    //indices from 1 to N-1 as we only want the inner vorticities
    double * v_test=new double[nsv];
    //print v_inner to check 
    cblas_dcopy(nsv,v_inner,1,v_test,1);
    
    
    for (int i=0;i<Ny-2;++i) {
        for (int j=0;j<Nx-2;++j){
        //taken as normal array now
            cout<<setw(4)<<v_test[i+Ny*j]<<"  ";
            
        }
        cout<<endl;
    }

    cout<<endl;

    for (int i=0;i<(Ny-2)*(Nx-2);++i) {
        cout<<setw(4)<<v_test[i]<<"  ";
    }
    cout<<endl;

     
   
   
    //for loop to check if correct
    for (int i=1;i<Nx-1;++i) {//columns 
        for (int j=1;j<Ny-1;++j){//rows
            v[i*Ny+j]=((one_dx2)*(s[(i+1)*Ny+j]-2*s[i*Ny+j]+s[(i-1)*Ny+j])+
            ((one_dy2)*(s[i*Ny+j+1]-2*s[i*Ny+j]+s[i*Ny+j-1])));
        }
    }
    cout<<endl;
    
     //print v again to check for loop solution
    for (int i=0;i<Ny;++i) {
        for (int j=0;j<Nx;++j){
        //taken as normal array now
            cout<<setw(4)<<v[i+Ny*j]<<"  ";   
        }
        cout<<endl;
    }
     cout<<endl;

       for (int i=0;i<(Ny)*(Nx);++i) {
        cout<<setw(4)<<v[i]<<"  ";
    }
    cout<<endl;
    //commplete vorticity matrix with boundaries 
delete [] v,s,v_inner,s_inner,A,v_test;
}
/*
 for (int i=0;i<ldA;++i) {
        for (int j=0;j<nsv;++j){
        //taken as normal array now
            cout<<setw(4)<<A[i+ldA*j]<<"  ";
            
        }
        cout<<endl;
    }
*/