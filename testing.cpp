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
    int Nx=10;//columns
    int Ny=10;//rows
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
    
    double *s_sub=s+Ny+1; 
    /*
    for (int i=0;i<Nx;++i) {
      s[(Ny*i)+Ny-2]=2.0;
    }
    */
    /*
    boundary consitions will be only used in Poisson solver RHS 
    keep convention of handout i.e i=columns, j=rows 
    */
    //vorticity. boundary conditions @ corners always 0..
    double* v=new double[(Nx)*(Ny)];
    /*top BC
    needs access to N_y-2 row of s-> rows*index + no. positions to last element
    */
    double* v_top=new double[Nx];
    /*bottom BC
    now we need row no. 1 on top (2nd row) of s -> rows*index+1 
    */
    double* v_bottom=new double[Nx];

    for (int i=0;i<Nx;++i){//Nx no of columns 
        v[Ny*i+Ny-1]=-s[(Ny*i)+Ny-2]*(2/pow(dy,2))-2*U/dy;//last row of v 
        v_top[i]=-s[(Ny*i)+Ny-2]*2/pow(dy,2)-2*U/dy;
        
        v[(Ny*i)]=-s[(Ny*i)+1]*2/pow(dy,2);//first row of v
        v_bottom[i]=-s[(Ny*i)+1]*2/pow(dy,2);
    }
    
    /**
    left needs access to 2nd (left) column 
    right to last to final column (Nx-2) s-> rows*columns 
    */
    double* v_left=new double[Ny];
    double* v_right=new double[Ny];


    for (int i=0;i<Ny;++i){
        v[i]=-s[Ny+i]*2/pow(dx,2);//first column of v 
        v_left[i]=-s[Ny+i]*2/pow(dx,2);
        v[(Nx-1)*Ny+i]=-s[(Nx-2)*Ny+i]*2/pow(dx,2);//last column of v
        v_right[i]=-s[(Nx-2)*Ny+i]*2/pow(dx,2);
    }


    //print s to check boundary conditions 
      for (int i=0;i<Ny;++i) {
        for (int j=0;j<Nx;++j){
        //taken as normal array now
            cout<<setw(4)<<s[i+Ny*j]<<"  ";
            
        }
        cout<<endl;
    }

    cout<<endl;

     double *s_inner=new double [(Ny-2)*(Nx-2)];
     for (int i=1;i<Nx-1;++i) {
        for (int j=1;j<Ny-1;++j){
        //taken as normal array now
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
    

    //print v
    for (int i=0;i<Ny;++i) {
        for (int j=0;j<Nx;++j){
        //taken as normal array now
            cout<<setw(4)<<v[i+Ny*j]<<"  ";
            
        }
        cout<<endl;
    }

    cout<<endl;

    delete [] v,s,s_sub,s_inner;
}
/* original A matrix
   A[0*ldA+3]=two_dxdy;
    A[1*ldA+3]=two_dxdy;
    A[2*ldA+3]=two_dxdy;
    //1st KU
     A[1*ldA+2]=one_dy2;
     A[2*ldA+2]=one_dy2;
     //2nd Ku
     A[2*ldA+1]=0.0;
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

    print_matrix(A,ldA,nsv);
*/