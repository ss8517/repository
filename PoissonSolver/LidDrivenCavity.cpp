#include<iostream>
#include<iomanip>
#include <cstdlib>
#include <fstream>
#include <string>   
#include <math.h>
#include <cblas.h>
#include <cmath>
#include "cstring"
#include<fstream>
#include "mpi.h"

#include "LidDrivenCavity.h"
#include "PoissonSolver.h"

using namespace std;

extern "C" {
    void Cblacs_pinfo(int*, int*);
    void Cblacs_get(int, int, int*);
    void Cblacs_gridinit(int*, char const*, int, int);
    void Cblacs_gridinfo(int, int*, int*, int*, int*);
    void Cblacs_barrier(int , char*);
    void Cblacs_gridexit(int);
    void Cblacs_exit(int);
}



LidDrivenCavity::LidDrivenCavity(int size,int nx,int ny) 
: size_(size),Nx(nx),Ny(ny)
{   
    s_inner=new double [size_];
        for (int i=0;i<(Nx-2)*(Ny-2);++i) {
        s_inner[i]=1.0;
    }
    v_inner=new double[size_]; 
    v_BC_top=new double[size_];
    v_BC_bottom=new double[size_];
    v_BC_left=new double[size_];
    v_BC_right=new double[size_];
    v_inner_col=new double[size_];
    s_inner_col=new double[size_];
    s_inner_row=new double[size_];
    s_inner_row_result1=new double[size_];
    s_inner_row_result2=new double[size_];
    v_inner_row=new double[size_];
    v_inner_row_result1=new double[size_];
    v_inner_row_result2=new double[size_];
    mult_1=new double[size_];
    mult_2=new double[size_];
    mult_3=new double[size_]; 
    

    ldA=3+(Ny-4);
    ldB=3;
    nsv=(Ny-2)*(Nx-2);
    
    A=new double [ldA*nsv];
    B_col=new double[ldB*nsv];
    B_row=new double[ldB*nsv];
    C_sb=new double [ldA*nsv];
    
    double* A_loc;
    double* B_loc;
    double*B_loc_red=new double[size_];

    int*    ipiv;
    int LWORK;
    double* WORK;
    double* AF;

    int* desca;
    int* descb;

    PoissonSolver* psolver=new PoissonSolver();

}

LidDrivenCavity::~LidDrivenCavity()
{   
    delete [] s_inner;
    delete [] v_inner;
    delete [] v_BC_top,v_BC_bottom,v_BC_left,v_BC_right;
    delete [] A,B_col,B_row,C_sb;
    delete [] s_inner_col;
    delete [] v_inner_col;
    delete [] v_inner_row,v_inner_row_result1,v_inner_row_result2;
    delete [] s_inner_row,s_inner_row_result1,s_inner_row_result2;
    delete [] mult_1,mult_2,mult_3;
    delete [] B_loc_red,A_loc,B_loc;
    delete [] desca,descb;
    delete [] ipiv,WORK,AF;
    delete psolver;
}

void LidDrivenCavity::SetDomainSize(double xlen, double ylen)
{
    Lx=xlen;
    Ly=ylen;
    cout << "lx is" <<Lx <<endl;
}

void LidDrivenCavity::SetGridSize(int nx, int ny)
{
    Nx=nx;
    Ny=Ny;
}

void LidDrivenCavity::SetTimeStep(double deltat)
{
    dt=deltat;
}

void LidDrivenCavity::SetFinalTime(double finalt)
{
    T=finalt;
}

void LidDrivenCavity::SetReynoldsNumber(double re)
{
    Re=re;
}

void LidDrivenCavity::SetPartitionSize(int px, int py)
{
    Px = px;
    Py = py;
}

void LidDrivenCavity::SetConstants () {
    //size_=(Ny-2)*(Nx-2);
    dx=Lx/(Nx-1.0);
    dy=Ly/(Ny-1.0);
    KL=2+(Ny-4);
    KU=KL;
    one_dx2=-1/dx/dx;
    one_dy2=-1/dy/dy;
    two_dxdy=2/dx/dx+2/dy/dy;

    one_2dy=1.0/2.0/dy;
    one_2dx=1.0/2.0/dx;

    //just for poisson
    ldAgb=1+2*KL+2*KU;//parallel
    one_Redx2=1/Re/dx/dx;
    one_Redy2=1/Re/dy/dy;
    two_Redxdy=-2/Re/dx/dx-2/Re/dy/dy;

}




void LidDrivenCavity::BoundaryConditions(){
    
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


void LidDrivenCavity::FillSBMatrix(int nsv, double* A, int ldA,int Nx,int Ny, double one_dx2,
double one_dy2,double two_dxdy){
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
void LidDrivenCavity::FillGBmatrix_1stDerCentral_col(){  
    for (int i=0; i<nsv; i++){
        if (i%(Ny-2)==0) {//every 3rd position 0 
            B_col[i*ldB]   =0.0;
        }
        else {
           B_col[i*ldB]   =one_2dy;
        }
        B_col[i*ldB+1] = 0.0; 
        if ((i+1)%(Ny-2)==0) {//every 3rd position 0 
            B_col[i*ldB+2]   =0.0;
        }
        else {
           B_col[i*ldB+2] = -one_2dy;
        }
    }
}

void LidDrivenCavity::FillGBmatrix_1stDerCentral_row(){  
    for (int i=0; i<nsv; i++){
        if (i%(Nx-2)==0) {//every 3rd position 0 
            B_row[i*ldB]   =0.0;
        }
        else {
           B_row[i*ldB]   =one_2dx;
        }
        B_row[i*ldB+1] = 0.0; 
        if ((i+1)%(Nx-2)==0) {//every 3rd position 0 
            B_row[i*ldB+2]   =0.0;
        }
        else {
           B_row[i*ldB+2] = -one_2dx;
        }
    }
}

void LidDrivenCavity::print_matrix (const double *A,const int row, const int col) {
    for (int i=0;i<row;++i) {
        for (int j=0;j<col;++j){
            cout<<setw(7)<<A[i+row*j]<<"  ";   
        }
        cout<<endl;
    }
    cout<<endl;
}

void LidDrivenCavity::print_matrix_row (const double *A,const int row, const int col) {//print matrix that is in row major format
    for (int i=0;i<row;++i) {
        for (int j=0;j<col;++j){
        //taken as normal array now
            cout<<setw(7)<<A[i*col+j]<<"  ";   
        }
        cout<<endl;
    }
    cout<<endl;
}


void LidDrivenCavity::InnerVorticity() {
     cblas_dsbmv(CblasColMajor, CblasUpper,nsv,KL,1.0,A,ldA,
                    s_inner,1,0.0,v_inner,1);
}

void LidDrivenCavity::Col2Row(double* A_col,double* B_row,int row,int col) {
    for (int i=0;i<col;++i) {
        for(int j=0;j<row;++j){
            B_row[col*j+i]=A_col[row*i+j];
        }
    }
}

void LidDrivenCavity::Row2Col(double* A_row,double* B_col,int row,int col) {
    for (int i=0;i<col;++i) {
        for(int j=0;j<row;++j){
            B_col[row*i+j]=A_row[col*j+i];
        }
    }
}


void LidDrivenCavity::Initialise(double xlen, double ylen,int nx, int ny,int px, int py, double deltat, double finalt, double re)
{
    this->SetDomainSize(xlen,ylen);
    this->SetGridSize(nx,ny);
    this->SetPartitionSize(px,py);
    this->SetTimeStep(deltat);
    this->SetFinalTime(finalt);
    this->SetReynoldsNumber(re);
    this->SetConstants();
    //this->SetArrays();
    //Initialise the Boundary Conditions
    this->BoundaryConditions();
    cout<<"The Bopundary conditions are: \n";
    cout<<"Top: \n";
    this->print_matrix(v_BC_top,Ny-2,Nx-2);
    cout<<"Bottom: \n";
    this->print_matrix(v_BC_bottom,Ny-2,Nx-2);
    cout<<"Left: \n";
    this->print_matrix(v_BC_left,Ny-2,Nx-2);
    cout<<"Right: \n";
    this->print_matrix(v_BC_right,Ny-2,Nx-2);
    //Fill Symmetric Banded matrix for inner vorticity calc
    this->FillSBMatrix(nsv,A,ldA,Nx,Ny,one_dx2,one_dy2,two_dxdy);
    //Inner Vorticity Calculation
    this->InnerVorticity();
    cout<<"Inner vorticity @ t \n";
    this->print_matrix(v_inner,Ny-2,Nx-2);

    //Fill the matrices for first two multiplications 
    //of first derivative central difference scheme
    //for vorticity @ t + dt
    
    this->FillGBmatrix_1stDerCentral_col();
    this->FillGBmatrix_1stDerCentral_row();

    cout<<"1st mult s_inner_col: \n";
    cblas_dgbmv(CblasColMajor,CblasNoTrans,nsv,nsv,1,1,1.0,B_col,ldB,s_inner,1,0.0,s_inner_col,1);
    print_matrix(s_inner_col,Ny-2,Nx-2);

    cout<<"1st mult v_inner_row_result2: \n";
    this->Col2Row(v_inner,v_inner_row,Nx-2,Ny-2);
    cblas_dgbmv(CblasColMajor,CblasNoTrans,nsv,nsv,1,1,1.0,B_row,ldB,v_inner_row,1,0.0,v_inner_row_result1,1);
    this->Row2Col(v_inner_row_result1,v_inner_row_result2,Ny-2,Nx-2);//v_inner_row_result2 is in column form
    this->print_matrix(v_inner_row_result2,Ny-2,Nx-2);
    //Add the boundary conditions
    //top and right added, left and bottom subtracted
    cout<<"1st mult V_inner_row_result2 after BCs added: \n";
    cblas_daxpy((Nx-2)*(Ny-2),-one_2dx,v_BC_left,1,v_inner_row_result2,1);//v_inner updated
    cblas_daxpy((Nx-2)*(Ny-2), one_2dx,v_BC_right,1,v_inner_row_result2,1);//
    this->print_matrix(v_inner_row_result2,Ny-2,Nx-2);
    

    cout<<"Mult_1: s_inner_col * v_inner_row_result2 : \n";
    cblas_dsbmv(CblasColMajor, CblasUpper,nsv, 0, 1.0, s_inner_col, 1,v_inner_row_result2,1, 0.0,mult_1, 1);
    this->print_matrix(mult_1,Ny-2,Nx-2);

    cout<<"2nd mult s_inner_row_result2: \n";
    this->Col2Row(s_inner,s_inner_row,Nx-2,Ny-2);
    cblas_dgbmv(CblasColMajor,CblasNoTrans,nsv,nsv,1,1,1.0,B_row,ldB,s_inner_row,1,0.0,s_inner_row_result1,1);
    this->Row2Col(s_inner_row_result1,s_inner_row_result2,Ny-2,Nx-2);
    this->print_matrix(s_inner_row_result2,Ny-2,Nx-2);
    
    cout<<"2nd mult v_inner_col: \n";
    cblas_dgbmv(CblasColMajor,CblasNoTrans,nsv,nsv,1,1,1.0,B_col,ldB,v_inner,1,0.0,v_inner_col,1);
    this->print_matrix(v_inner_col,Ny-2,Nx-2);

    cout<<"2nd mult v_inner_col  after BCs added: \n";
    cblas_daxpy((Nx-2)*(Ny-2),-one_2dy,v_BC_bottom,1,v_inner_col,1);
    cblas_daxpy((Nx-2)*(Ny-2),one_2dy,v_BC_top,1,v_inner_col,1);
    this->print_matrix(v_inner_col,Ny-2,Nx-2);

    cout<<"Mult_2: s_inner_row_result2 * v_inner_col: \n";
    cblas_dsbmv(CblasColMajor, CblasUpper,nsv, 0, 1.0, s_inner_row_result2, 1,v_inner_col,1, 0.0,mult_2, 1);
    this->print_matrix(mult_2,Ny-2,Nx-2);

    //Multiplication 3
    //second der central difference scheme 
    cout<<"Mult_3: C_sb * v_inner \n";
    this->FillSBMatrix(nsv, C_sb,ldA,Nx,Ny,one_Redx2,one_Redy2,two_Redxdy);
    cblas_dsbmv(CblasColMajor, CblasUpper,nsv,KL,1.0,C_sb,ldA,
                    v_inner,1,0.0,mult_3,1);
    //Boundary conditions Added 
    cblas_daxpy((Nx-2)*(Ny-2),one_Redy2,v_BC_bottom,1,mult_3,1);
    cblas_daxpy((Nx-2)*(Ny-2),one_Redy2,v_BC_top,1,mult_3,1);
    cblas_daxpy((Nx-2)*(Ny-2),one_Redx2,v_BC_left,1,mult_3,1);
    cblas_daxpy((Nx-2)*(Ny-2),one_Redx2,v_BC_right,1,mult_3,1);
    this->print_matrix(mult_3,Ny-2,Nx-2);

    
    //Adding the results of the 3 multiplications 
    //to get vorticity @ t + dt
    cout<<"Addition Result:  \n";
    cblas_daxpy((Ny-2)*(Nx-2),-1.0,mult_1,1,mult_3,1);
    cblas_daxpy((Ny-2)*(Nx-2),1.0,mult_2,1,mult_3,1);
    print_matrix(mult_3,Ny-2,Nx-2);

    cout<<"v_inner @ t + Dt \n";
    cblas_dscal((Ny-2)*(Nx-2),dt,mult_3,1);
    cblas_daxpy((Ny-2)*(Nx-2),1.0,v_inner,1,mult_3,1);
    print_matrix(mult_3,Ny-2,Nx-2);


    psolver->BuildLocalMatrices( s_inner,A_loc,  B_loc,B_loc_red, mult_3, size_, ldAgb, Nx, Ny, one_dx2, one_dy2, two_dxdy,ipiv, KL, KU
        , desca, AF, WORK, descb) ;
    
    //print_matrix(A_loc,ldAgb,NB);

    //psolver->LUfactorisation(A_loc,N,ipiv,BWU,BWL,desca,AF,LAF,WORK,LWORK,mycol);
    
    //psolver->LinearSolver(A_loc,B_loc,N,ipiv,BWU,BWL,desca,descb,AF,LAF,WORK,LWORK,mycol);

 
 
}