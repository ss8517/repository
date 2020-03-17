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

void Col2Row(double* A_col,double* B_row,int row,int col) {
    for (int i=0;i<col;++i) {
        for(int j=0;j<row;++j){
            B_row[col*j+i]=A_col[row*i+j];
        }
    }
}

int main () {
    double U=1.0;
    double Lx=1.0;
    double Ly=1.0;
    double Re=100;
    int Nx=4;//columns
    int Ny=6;//rows
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
    


    double *s_inner=new double [(Ny-2)*(Nx-2)];
    double *s_inner_test=new double [(Ny-2)*(Nx-2)];
    double *s_inner_test_2=new double [(Ny-2)*(Nx-2)];
   for (int i=0;i<(Nx-2)*(Ny-2);++i) {
        s_inner[i]=i;
    }
    
   print_matrix(s_inner,(Ny-2),(Nx-2));
    Col2Row(s_inner,s_inner_test,Ny-2,Nx-2);
    Row2Col(s_inner_test,s_inner_test_2,Ny-2,Nx-2);
    
    for (int i=0;i<(Nx-2)*(Ny-2);++i) {
        cout<< s_inner_test[i]<<" ";
    }
    cout<<endl;
    
    for (int i=0;i<(Nx-2)*(Ny-2);++i) {
        cout<< s_inner[i]<<" ";
    }
    cout<<endl;
    for (int i=0;i<(Nx-2)*(Ny-2);++i) {
        cout<< s_inner_test_2[i]<<" ";
    }

cout<<endl;
double* v_BC_top=new double[(Nx-2)*(Ny-2)];
memset(v_BC_top, 0.0, (Nx-2)*(Ny-2)*sizeof(v_BC_top[0]));

double* v_BC_bottom=new double[(Nx-2)*(Ny-2)];
memset(v_BC_bottom, 0.0, (Nx-2)*(Ny-2)*sizeof(v_BC_bottom[0]));

double* v_BC_left=new double[(Nx-2)*(Ny-2)];
memset(v_BC_left, 0.0, (Nx-2)*(Ny-2)*sizeof(v_BC_left[0]));

double* v_BC_right=new double[(Nx-2)*(Ny-2)];
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

 print_matrix(v_BC_top,Ny-2,Nx-2);
 print_matrix(v_BC_bottom,Ny-2,Nx-2);
 print_matrix(v_BC_left,Ny-2,Nx-2);
 print_matrix(v_BC_right,Ny-2,Nx-2);



 delete [] s,s_inner,v_BC_top,v_BC_bottom,v_BC_left,v_BC_right;

}


   
    /*
    for (int i=1;i<Nx-1;++i) {
        for (int j=1;j<Ny-1;++j){
        //taken as normal array now
           s_inner[(i-1)*Ny+j-1]=s[i*Ny+j];
        }
     }

    /*top BC
    needs access to N_y-2 row of s-> rows*index + no. positions to last element
    */

    /*
    boundary consitions will be only used in Poisson solver RHS 
    keep convention of handout i.e i=columns, j=rows 
    */
    //vorticity. boundary conditions @ corners v always 0
    
     
    
    /*bottom BC
    now we need 2nd row of s -> rows*index+1 
     for (int i=0;i<Ny;++i){
        //v_BC[i]=-s[Ny+i]*2/dx/dx;//first column of v (left)
        //v_BC[(Nx-1)*Ny+i]=-s[(Nx-2)*Ny+i]*2/dx/dx;//last column of v (right)
        v[i]=-s[Ny+i]*2/dx/dx;//first column of v (left)
        v[(Nx-1)*Ny+i]=-s[(Nx-2)*Ny+i]*2/dx/dx;//last column of v (right)

    }
    */
            
/*
    Boundaries
      for (int i=0;i<Nx;++i){//Nx no of columns 
        //v_BC[Ny*i+Ny-1]=-s[(Ny*i)+Ny-2]*2/dy/dy-2*U/dy;//last row of s
        v[Ny*i+Ny-1]=-s[(Ny*i)+Ny-2]*2/dy/dy-2*U/dy;//last row of s
        //v_BC[(Ny*i)]=-s[(Ny*i)+1]*2/dy/dy;//first row of v
        v[(Ny*i)]=-s[(Ny*i)+1]*2/dy/dy;//first row of v
    }

  //for loops for check
    for (int i=1;i<Nx-1;++i) {
        for (int j=1;j<Ny-1;++j){//rows
            v[i*Ny+j]=((one_dx2)*(s[(i+1)*Ny+j]-2*s[i*Ny+j]+s[(i-1)*Ny+j])+
            ((one_dy2)*(s[i*Ny+j+1]-2*s[i*Ny+j]+s[i*Ny+j-1])));
        }
    }
    cout<<"Inner vorticity @ t (for loops):  \n";
    print_matrix(v,Ny,Nx);

  //stream function initialisation 
   
    double* s=new double[(Nx)*(Ny)]; 
    memset(s, 0, (Nx)*(Ny)*sizeof(s[0]));//to fill array with 0s - check it 
    //initialise stream function once

    for (int i=1;i<Nx-1;++i) {
        for (int j=1;j<Ny-1;++j){
            s[i*Ny+j]=1.0;
        }
    }
    print_matrix(s,Ny,Nx);
 for (int i=0;i<(Nx-1)*(Ny-1);++i) {
        s_inner[i]=1;
    }
     */
   