#include "LidDrivenCavity.h"
#include "PoissonSolver.h"
#include <iostream>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <string>
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

LidDrivenCavity::LidDrivenCavity()
{
}

//de-structor - deallocate memory
LidDrivenCavity::~LidDrivenCavity()
{
    delete[] psi;
    delete[] omega;
    delete[] A;
    delete[] psi_in;
    delete[] omega_in;
    delete[] A_Nj;
    delete[] A_Ni;
    delete[] res_psi_j;
    delete[] res_psi_i;
    delete[] res_omega_j;
    delete[] res_omega_i;
    delete[] res_omega_in3;
    delete[] omega_in_new;
    delete[] mult1;
    delete[] mult2;
    delete[] omega_in_top;
    delete[] omega_in_bottom;
    delete[] omega_in_left;
    delete[] omega_in_right;
    delete[] psi_new;
    delete[] A_generalbanded;
    delete[] ipiv;
    delete[] psi_new;
}

void LidDrivenCavity::SetDomainSize(double xlen, double ylen)
{
    Lx = xlen;
    Ly = ylen;
}

void LidDrivenCavity::SetGridSize(int nx, int ny)
{
    Nx = nx;
    Ny = ny;
}

void LidDrivenCavity::SetPartitionSize(int px, int py)
{
    Px = px;
    Py = py;
}

void LidDrivenCavity::SetTimeStep(double deltat)
{
   dt = deltat;
}

void LidDrivenCavity::SetFinalTime(double finalt)
{
   T = finalt;
}

void LidDrivenCavity::SetReynoldsNumber(double re){
    Re = re;
}

void LidDrivenCavity::SetCoeffsandsizes(){
     dx = Lx / (Nx - 1.0);
     dy = Ly / (Ny - 1.0);
     lda_Nj = 3;
     lda = 3+Ny-4; // for the first A matrix upper symmetric
     nsv = (Nx-2)*(Ny-2);
     lda2 = (1+2*(Ny-2)+(Ny-2));
     nsv2 = (Ny-2)*(Nx-2);
     U = 1.0;
     
 }
 
void LidDrivenCavity:: SetBCdata(){
    omega_in_top = new double[nsv];
    omega_in_bottom = new double[nsv];
    omega_in_left = new double[nsv];
    omega_in_right = new double[nsv];
}

void LidDrivenCavity::SetConstMatrices(){
    //parameters for omega at t+dt
    A = new double[lda*nsv]; // A matrix upper symmetric  initialisation
    A_generalbanded = new double[lda2*nsv2];
    A_Nj = new double[lda_Nj*nsv];
    A_Ni = new double[lda_Nj*nsv];
    
}

void LidDrivenCavity::SetVectors(){
    
    psi = new double[Nx*Ny];
    omega = new double[Nx*Ny];
    omega_in = new double[nsv]; //interior points vectors
    psi_in = new double[nsv];
    res_psi_j = new double[nsv];         //column result
    res_psi_i = new double[nsv];         //row result
    res_omega_j = new double[nsv];       //row result
    res_omega_i = new double[nsv];       //column result
    res_omega_in3 = new double[nsv];      //result of the third multiplication
    omega_in_new = new double[nsv]; 
     
   //multiplication results
    mult1 = new double[nsv];
    mult2 = new double[nsv];
    
    psi_new = new double[nsv]; //result of the poisson
    //Initialise omega and psi to 0.0 at t = 0, initial conditions
    ipiv = new int[nsv2];
    
    for (int i = 0; i<Nx; i++){
        for (int j = 0; j<Ny; j++){
            psi[i*Ny+j] = 0.0;
            omega[i*Ny+j] = 0.0;
        }
    }
    
     for (int i=0;i<(Nx-1)*(Ny-1);i++) {
        psi_in[i]=1.0;
        omega_in[i]=1.0;
    }
    
     for (int i = 0; i<Nx-2; i++){
        for (int j = 0; j<Ny-2; j++){
            ipiv[i*(Ny-2)+j] = 0.0;
        }
    }
   
}

void LidDrivenCavity::PrintMatrix(int Ny, int Nx, double* H) {
for (int i = 0; i < Nx; ++i) {
    for (int j = 0; j < Ny; ++j) {
        cout << setw(5) << H[j*Nx+i] << " ";
    }
        cout << endl;
}  
}  

//Vorticity boundary conditions implementation at time t.
//Top boundary conditions:
void LidDrivenCavity::TopBC(){
    for (int i=0; i<Nx; i++){
        omega[Ny*i+Ny-1] = (-psi[(Ny*i)+Ny-2]*2.0/(pow(dy,2.0))-2*U/dy);
    }
}

//Bottom boundary conditions
void LidDrivenCavity::BottomBC(){
    for (int i=0; i<Nx; i++){
      omega[(Ny*i)] = (-psi[(Ny*i)+1])*2.0/(pow(dy,2.0));
    }
}

//Left boundary conditions
void LidDrivenCavity::LeftBC(){
    for (int j=0; j<Ny; j++){
        omega[j] = (-psi[j+Ny])*2.0/(pow(dx,2.0));
    }
}

//Right boundary conditions
void LidDrivenCavity::RightBC(){
    for (int j=0; j<Ny; j++){
        omega[(Nx-1)*Ny+j] = (-psi[(Nx-2)*Ny+j])*2.0/(pow(dx,2.0));
   }
}

//______________________________________________________________________________
void LidDrivenCavity::AddBC(){
    
    //New vectors for each BC that needs to be added to the central differrences
    omega_in_top = new double[nsv];
    omega_in_bottom = new double[nsv];
    omega_in_left = new double[nsv];
    omega_in_right = new double[nsv];
    
    //initialise
    for (int i=0;i<Nx-2;i++){
        for (int j=0; j<Ny-2; j++){
            omega_in_top[i*(Ny-2)+j]=0.0;
            omega_in_bottom[i*(Ny-2)+j]=0.0;
            omega_in_left[i*(Ny-2)+j]=0.0;
            omega_in_right[i*(Ny-2)+j]=0.0;
        }
    }
   
    for (int i=0;i<Nx-2;i++){
        omega_in_top[(Ny-2)*i+Ny-3] = -psi_in[(Ny-2)*i+Ny-3]*2/pow(dy,2)-2*U/dy;
        omega_in_bottom[(Ny-2)*i] = -psi_in[(Ny-2)*i]*2/pow(dy,2);
    }
    
    for (int j=0;j<Ny-2;j++){
        omega_in_left[j] = -psi_in[j]*2/(pow(dx,2));
        omega_in_right[(Nx-3)*(Ny-2)+j] = -psi_in[(Nx-3)*(Ny-2)+j]*2/(pow(dx,2));
    }
    
}

//______________________________________________________________________________
//Transform from Column to Row Major function
void LidDrivenCavity::TransColtoRow(double *A, int Ny, int Nx){
    double *temp = new double[Ny*Nx];
    for (int i=0; i<Nx; i++){
        for (int j = 0; j<Ny; j++){
            temp[j*Nx+i] = A[i*Ny+j];
            //B[Nx*j+i] = A[i*Ny+j];
            }
    }
    cblas_dcopy(Nx*Ny,temp,1,A,1);
    delete []temp;
}

//Transform from Row to Column Major form
void LidDrivenCavity::TransRowtoCol(double *A, int Ny, int Nx){
    double *temp = new double[Ny*Nx];
    for (int i=0; i<Nx; i++){
        for (int j = 0; j<Ny; j++){
            temp[i*Ny+j] = A[j*Nx+i];
        }
    }
    cblas_dcopy(Nx*Ny,temp,1,A,1);
    delete []temp;
}

//_______________________________________________________________________________
//Create Banded Upper symmetric matrix A for omega calculation at time t
void LidDrivenCavity::Amatrix(){
    double coeffx = -1.0/(pow(dx,2.0));
    double coeffy = -1.0/(pow(dy,2.0));
    double diag = 2.0/(pow(dx,2.0))+ 2.0/(pow(dy,2.0));
    int nsv = (Nx-2)*(Ny-2);
    int lda = 3+Ny-4;
     for (int i = 0; i < nsv; i++) {
          A[i*lda]=coeffx;
          for (int j=1;j<(Ny-3);++j){
            A[i*lda+j]=0.0;
          }
          if (i==0 || i%(Ny-2)==0) {
              A[i*lda+(Ny-4)+1]=0.0;
            }
          else {
               A[i*lda+(Ny-4)+1]=coeffy;
          }
          A[i*lda+(Ny-4)+2]=diag;
        }
        PrintMatrix(Ny-2,Nx-2,A);
}

//Create general Banded Upper symmetric matrix A
void LidDrivenCavity::AmatrixGB(){
    double coeffx = -1.0/(pow(dx,2.0));
    double coeffy = -1.0/(pow(dy,2.0));
    double diag = 2.0/(pow(dx,2.0))+ 2.0/(pow(dy,2.0));      
   
    for (int i = 0; i < nsv; ++i) {
        A_generalbanded[i*lda2+(Ny-2)]=coeffx;//doesnt change - final Kl row
        for (int j=1;j<(Ny-3);++j) {//populate all diagonals that have 0s
            A_generalbanded[i*lda2+j+(Ny-2)]=0.0;
        }
        if (i%(Ny-2)==0) {//every 3rd position 0
            A_generalbanded[i*lda2+(Ny-4)+1+(Ny-2)]=0.0;
        }
        else {
            A_generalbanded[i*lda2+(Ny-4)+1+(Ny-2)]=coeffy;
        }
        A_generalbanded[i*lda2+(Ny-4)+2+(Ny-2)]=diag;//row of 4s
        if ((i+1)%(Ny-2)==0) {//every third position 0, offset 1
            A_generalbanded[i*lda2+(Ny-4)+3+(Ny-2)]=0.0;
        }
        else {
            A_generalbanded[i*lda2+(Ny-4)+3+(Ny-2)]=coeffy;
        }
        for (int j=1;j<(Ny-3);++j) {//populate all diagonals that have 0s
            A_generalbanded[i*lda2+(Ny-4)+3+j+(Ny-2)]=0.0;
        }
        A_generalbanded[i*lda2+(2*(Ny-4)+3)+1+(Ny-2)]=coeffx;
    }
    cout <<"A_generalbanded" <<endl;
    PrintMatrix(Ny-2,Nx-2,A_generalbanded);
}
   
/* Interior vorticity omega at t+dt:
Only accounting interior points initially, then adding BCS
psi_in, psi_omega are initially considered
A_Nj - A matrix of dimensions (Ny-2)*(Nx-2)
A_Ni - A matrix of dimensions (Nx-2)*(Ny-2)
*/    
void LidDrivenCavity::A_Nj_matrix(){
    double coeffy = 1.0/(2.0*dy);
//where j stands for column
// add zero every third positions in the Kl and Ku

    for (int i=0; i<nsv; i++){
        if (i%(Ny-2)==0){
            A_Nj[i*lda_Nj] = 0.0;
        }
        else{
            A_Nj[i*lda_Nj] = coeffy;
        }
        A_Nj[i*lda_Nj+1] = 0.0;
        if ((i+1)%(Ny-2)==0) {
            A_Nj[i*lda_Nj+2] = 0.0;
        }
        else{
            A_Nj[i*lda_Nj+2] = -coeffy;
        }
    }
    
    cout <<"A_Nj_matrix" <<endl;
    PrintMatrix(Ny-2,Nx-2,A_Nj);
}

void LidDrivenCavity::A_Ni_matrix(){
    double coeffx = 1.0/(2.0*dx);

    for (int i=0; i<nsv; i++){
        if (i%(Nx-2)==0){
            A_Ni[i*lda_Nj] = 0.0;
        }
        else{
            A_Ni[i*lda_Nj] = coeffx;
        }
        A_Ni[i*lda_Nj+1] = 0.0;
        if ((i+1)%(Nx-2)==0) {
            A_Ni[i*lda_Nj+2] = 0.0;
        }
        else{
            A_Ni[i*lda_Nj+2] = -coeffx;
        }
    }
}

//Function that multiplies the vectors of psi with omega
void LidDrivenCavity::Afirstmult(){
   //First multiplication result: A_Nj* psi_in:
   cout << "nsv: " << nsv <<endl;
   cout << "lda_Nj" <<lda_Nj << endl;
   cout <<"psi_in" << endl;
   PrintMatrix(Ny-2,Nx-2,psi_in);
       cout <<"A_Nj_matrix" <<endl;
    PrintMatrix(Ny-2,Nx-2,A_Nj);
    
    cblas_dgbmv(CblasColMajor,CblasNoTrans,nsv,nsv,1,1,1.0,A_Nj,lda_Nj,psi_in,1,0.0,res_psi_j,1);
    
    TransColtoRow(omega_in,Ny-2,Nx-2);
    cblas_dgbmv(CblasColMajor,CblasNoTrans,nsv,nsv,1,1,1.0,A_Ni,3,omega_in,1,0.0,res_omega_i,1); //Result before BCs
    TransRowtoCol(res_omega_i,Ny-2,Nx-2);

    
    //transform back omega inner.
    TransRowtoCol(omega_in,Ny-2,Nx-2);
    
    //Add Boundary condition vectors
    cblas_daxpy(nsv,-1/(2*dx),omega_in_left,1,res_omega_i,1); //left
    cblas_daxpy(nsv,1/(2*dx),omega_in_right,1,res_omega_i,1); //right
   
    //1st multiplication
    cblas_dsbmv(CblasColMajor,CblasUpper,nsv,0,1.0,res_psi_j,1,res_omega_i,1,0.0,mult1,1);    
    /*
    cout <<"Nx" << endl;
    cout << Nx << endl;
    cout << "A_Nj" <<endl;
    PrintMatrix(nsv,lda_Nj,res_omega_i);
    cout << "dx" << endl;
    cout << dx << endl;
    cout << "Res omega: "<< endl;
    PrintMatrix(Ny-2,Nx-2,res_omega_i);
     * */
   
}  

void LidDrivenCavity::Asecondmult(){
   
    TransColtoRow(psi_in,Ny-2,Nx-2);
    cblas_dgbmv(CblasColMajor,CblasNoTrans,nsv,nsv,1,1,1.0,A_Ni,lda_Nj,psi_in,1,0.0,res_psi_i,1);
    TransRowtoCol(res_psi_i,Ny-2,Nx-2);
    
    //transform back psi_in:
    TransRowtoCol(psi_in,Ny-2,Nx-2);
    
    cblas_dgbmv(CblasColMajor,CblasNoTrans,nsv,nsv,1,1,1.0,A_Nj,lda_Nj,omega_in,1,0.0,res_omega_j,1);

    cblas_daxpy(nsv,1/(2*dy),omega_in_top,1,res_omega_j,1); //top
    cblas_daxpy(nsv,-1/(2*dy),omega_in_bottom,1,res_omega_j,1);  //bottom
    
    cblas_dsbmv(CblasColMajor,CblasUpper,nsv,0,1.0,res_psi_i,1,res_omega_j,1,0.0,mult2,1);    
}

void LidDrivenCavity::Athirdmult(){

    cblas_dsbmv(CblasColMajor,CblasUpper,nsv,Ny-2,1.0,A,lda,omega_in,1.0,0.0,res_omega_in3,1.0);
    //omegat3(A,omega,res_omega_in3,psi,psi_in,dx,dy,lda,nsv,Ny,Nx);
    //add the Bcs    
    cblas_daxpy(nsv,1/(dx*dx),omega_in_left,1,res_omega_in3,1); //left
    cblas_daxpy(nsv,1/(dx*dx),omega_in_right,1,res_omega_in3,1); //right
    cblas_daxpy(nsv,1/(dy*dy),omega_in_top,1,res_omega_in3,1); //top
    cblas_daxpy(nsv,1/(dy*dy),omega_in_bottom,1,res_omega_in3,1);  //bottom
    cblas_dscal(nsv,1/Re,res_omega_in3,1);

}

void LidDrivenCavity::Finalcalc(){
    cblas_dcopy(nsv,mult1,1,res_omega_in3,1);
    cblas_daxpy(nsv,-1.0,mult2,1,mult1,1); //result saved in mult1
    
    cblas_daxpy(nsv,-1.0,mult1,1,res_omega_in3,1); //result saved in res_omega_in3
    cblas_dcopy(nsv,omega_in,1,omega_in_new,1);
    cblas_daxpy(nsv,dt,res_omega_in3,1,omega_in_new,1);
    cout << "New vorticity omega: " << endl;
    PrintMatrix(Ny-2,Nx-2,omega_in_new);
}

void LidDrivenCavity::omegat3(){
    int newlda = Ny-2;
    cblas_dsbmv(CblasColMajor,CblasUpper, nsv, newlda, 1.0,A, lda, psi_in, 1.0, 0.0,omega_in, 1.0);
    cout <<"omega at t=0 " <<endl;
    PrintMatrix(Ny-2,Nx-2,omega_in);
}


void LidDrivenCavity::Initialise(double xlen, double ylen,int nx, int ny,int px, int py, double deltat, double finalt, double re)
{
    //Initialise and Set all arguments
    this->SetDomainSize(xlen,ylen);
    this->SetGridSize(nx,ny);
    this->SetPartitionSize(px,py);
    this->SetTimeStep(deltat);
    this->SetFinalTime(finalt);
    this->SetReynoldsNumber(Re);
    this->SetCoeffsandsizes();
    //Initialise arrays, fixed matrices, sizes, coefficients needed
    
    //this->SetBCdata();
    this->SetConstMatrices();
    this->SetVectors();
    this->SetCoeffsandsizes();
    this->Amatrix();
    this->AmatrixGB();
    this->A_Nj_matrix();
    this->A_Ni_matrix();
    
    //Add boundary conditions
    this->TopBC();
    this->BottomBC();
    this->LeftBC();
    this->RightBC();
    
    this->AddBC();
    //compute omega at time t = 0
    this->omegat3();
    
    this->Afirstmult();
    this->Asecondmult();
    this->Athirdmult();
    this->Finalcalc();
    
 
    PoissonSolver* solver=new PoissonSolver();
     solver->LUfactorisation(double* A_generalbanded, int nsv2, int lda2, int* ipiv, int Ny, int Nx);
     solver->SolveLapack(double* A_generalbanded, double* omega_in_new, double *psi_new, int nsv2, int lda2, int* ipiv, int Ny, int Nx);
     

}


void LidDrivenCavity::Integrate()
{
    
     //compute initial omega t+dt 

     

}

