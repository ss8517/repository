#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <fstream>
using namespace std;
#include <string>   
#include <math.h>
#include <cblas.h>
#include <cmath>
#include "cstring"
#include "mpi.h"

#include "PoissonSolver.h"
#include "LidDrivenCavity.h"

#define F77NAME(x) x##_
extern "C" {
    
    void F77NAME(pdgbtrf) (const int& 	N,const int& 	BWL, const int& 	BWU, double *A, const int& JA,
    int* DESCA,int * IPIV, double *AF, const int& LAF, double  * WORK, const int&	LWORK,
    int * INFO);

    
    void F77NAME(pdgbtrs) (const char& TRANS, const int& N, const int& BWL, const int& BWU,
    const int& NRHS, double * A, const int& JA, int * DESCA,int * IPIV,
    double * B, const int& IB, int * DESCB, double * AF, const int& LAF,
    double * WORK, const int& LWORK, int * info);

    void Cblacs_pinfo(int*, int*);
    void Cblacs_get(int, int, int*);
    void Cblacs_gridinit(int*, char const*, int, int);
    void Cblacs_gridinfo(int, int*, int*, int*, int*);
    void Cblacs_barrier(int , char*);
    void Cblacs_gridexit(int);
    void Cblacs_exit(int);
}

PoissonSolver::PoissonSolver () {}

PoissonSolver::~PoissonSolver () {}

void PoissonSolver::print_matrix (const double *A,const int row, const int col) {
    for (int i=0;i<row;++i) {
        for (int j=0;j<col;++j){
            cout<<setw(7)<<A[i+row*j]<<"  ";   
        }
        cout<<endl;
    }
    cout<<endl;
}

void PoissonSolver::BuildLocalMatrices (double* s_inner,double* A_loc, double* B_loc,double* B_loc_red, double* mult_3,int size_,int ldAgb,int Nx,int Ny,double one_dx2,double one_dy2,double two_dxdy,int* ipiv,int KL,int KU,
        int* desca,double* AF,double* WORK,int* descb) 
{

    int myid,npe, nproc, ctxt, myrow, mycol;
    int nprow,npcol; //processes grid

    Cblacs_pinfo(&myid, &npe);
    Cblacs_get( -1, 0, &ctxt );
    Cblacs_gridinit( &ctxt,"Row_Major",1,npe);
    Cblacs_gridinfo( ctxt, &nprow, &npcol, &myrow, &mycol);
    
    int N=size_;
    int BWL=KL;
    int BWU=KU;
    int NB=ceil((double)N/npcol);
    int NRHS=1;
    int JA=1;
    int IB=1;
    int LA=(1 + 2*BWL + 2*BWU)*NB;//ldAgb*NB
    int LAF=(NB+BWU)*(BWL+BWU)+6*(BWL+BWU)*(BWL+2*BWU);
    //+max(NRHS*(NB+2*BWL+4*BWU), 1);
    ipiv = new int [N];//might have to be just int 
    int LWORK=LAF;
    WORK=new double [LWORK];
    AF=new double [LAF];
    
    
    desca=new int [7];
    descb=new int [7];

    desca[0] = 501;
    desca[1] = ctxt;
    desca[2] = N;
    desca[3] = NB;
    desca[4] = 0;
    desca[5] = 1+2*BWL+2*BWU;//local leading dim
    desca[6] = 0;

    

    descb[0] = 502;
    descb[1] = ctxt;
    descb[2] = N;
    descb[3] = NB;
    descb[4] = 0;
    descb[5] = NB;
    descb[6] = 0;


    for (unsigned int p=0;p<npcol;++p) {//npcol is the process number i.e column on cblacs grid 
        if (N%npcol!=0 && mycol==p && p==(npcol-1)) {
            cout<<"Last Process, nprow: "<<nprow<<" "<<"npcol: "<<npcol<<" "<<"myrow: "<<myrow<<" "<<"mycol: "<<mycol<<endl;
            
            A_loc=new double[ldAgb*(N-NB*(npcol-1))];
            B_loc=new double [(N-NB*(npcol-1))];//local matrix 
            
            for (int i = 0; i < (N-NB*(npcol-1)); ++i) {//size of col of last row is N-NB*(npcol-1) if N/p is not integer
                A_loc[i*ldAgb+2*(Ny-2)]=one_dx2;//doesnt change - final Kl row 
            
                for (int j=1;j<(Ny-3);++j) {//populate all diagonals that have 0s
                    A_loc[i*ldAgb+j+2*(Ny-2)]=0.0;
                }

                if ((i+p*NB)%(Ny-2)==0) {//every 3rd position 0 
                    A_loc[i*ldAgb+(Ny-4)+1+2*(Ny-2)]=0.0;
                }

                else {
                A_loc[i*ldAgb+(Ny-4)+1+2*(Ny-2)]=one_dy2;
                }

                A_loc[i*ldAgb+(Ny-4)+2+2*(Ny-2)]=two_dxdy;//row of 4s
        
                if (((i+1)+p*NB)%(Ny-2)==0) {//every third position 0, offset 1 
                    A_loc[i*ldAgb+(Ny-4)+3+2*(Ny-2)]=0.0;
                }
        
                else {
                    A_loc[i*ldAgb+(Ny-4)+3+2*(Ny-2)]=one_dy2;
                }

                for (int j=1;j<(Ny-3);++j) {//populate all diagonals that have 0s
                    A_loc[i*ldAgb+(Ny-4)+3+j+2*(Ny-2)]=0.0;
                }
                
                A_loc[i*ldAgb+(2*(Ny-4)+3)+1+2*(Ny-2)]=one_dx2;
                
                B_loc[i]=mult_3[i+NB*p];
            }
        }

        else if (N%npcol!=0 && mycol==p && p!=(npcol-1)) {
            cout<<"nprow: "<<nprow<<" "<<"npcol: "<<npcol<<" "<<"myrow: "<<myrow<<" "<<"mycol: "<<mycol<<endl;
                       
            A_loc =new double[LA];
            B_loc=new double [NB];//local matrix 
            
            for (int i = 0; i < NB; ++i) {
                A_loc[i*ldAgb+2*(Ny-2)]=one_dx2;//doesnt change - final Kl row 
            
                for (int j=1;j<(Ny-3);++j) {//populate all diagonals that have 0s
                    A_loc[i*ldAgb+j+2*(Ny-2)]=0.0;
                }

                if ((i+p*NB)%(Ny-2)==0) {//every 3rd position 0 //p*NB to divide matrix 
                    A_loc[i*ldAgb+(Ny-4)+1+2*(Ny-2)]=0.0;
                }

                else {
                    A_loc[i*ldAgb+(Ny-4)+1+2*(Ny-2)]=one_dy2;
                }

                A_loc[i*ldAgb+(Ny-4)+2+2*(Ny-2)]=two_dxdy;//row of 4s
        
                if (((i+1)+p*NB)%(Ny-2)==0) {//every third position 0, offset 1 
                    A_loc[i*ldAgb+(Ny-4)+3+2*(Ny-2)]=0.0;
                }
        
                else {
                    A_loc[i*ldAgb+(Ny-4)+3+2*(Ny-2)]=one_dy2;
                }

                for (int j=1;j<(Ny-3);++j) {//populate all diagonals that have 0s
                    A_loc[i*ldAgb+(Ny-4)+3+j+2*(Ny-2)]=0.0;
                }
        
                A_loc[i*ldAgb+(2*(Ny-4)+3)+1+2*(Ny-2)]=one_dx2;
                
                B_loc[i]=mult_3[p*NB+i];   
            }
        }
        
        else if (N%npcol==0 && mycol==p) {
            cout<<"nprow: "<<nprow<<" "<<"npcol: "<<npcol<<" "<<"myrow: "<<myrow<<" "<<"mycol: "<<mycol<<endl;
            
            A_loc =new double[LA];
            B_loc=new double [NB];//local matrix 
            
            for (int i = 0; i < NB; ++i) {
                
                A_loc[i*ldAgb+2*(Ny-2)]=one_dx2;//doesnt change - final Kl row 
            
                for (int j=1;j<(Ny-3);++j) {//populate all diagonals that have 0s
                    A_loc[i*ldAgb+j+2*(Ny-2)]=0.0;
                }

                if ((i+p*NB)%(Ny-2)==0) {//every 3rd position 0 //p*NB to divide matrix 
                    A_loc[i*ldAgb+(Ny-4)+1+2*(Ny-2)]=0.0;
                }

                else {
                    A_loc[i*ldAgb+(Ny-4)+1+2*(Ny-2)]=one_dy2;
                }

                A_loc[i*ldAgb+(Ny-4)+2+2*(Ny-2)]=two_dxdy;//row of 4s
        
                if (((i+1)+p*NB)%(Ny-2)==0) {//every third position 0, offset 1 
                    A_loc[i*ldAgb+(Ny-4)+3+2*(Ny-2)]=0.0;
                }
        
                else {
                    A_loc[i*ldAgb+(Ny-4)+3+2*(Ny-2)]=one_dy2;
                }

                for (int j=1;j<(Ny-3);++j) {//populate all diagonals that have 0s
                    A_loc[i*ldAgb+(Ny-4)+3+j+2*(Ny-2)]=0.0;
                }
        
                A_loc[i*ldAgb+(2*(Ny-4)+3)+1+2*(Ny-2)]=one_dx2;
        
                B_loc[i]=mult_3[p*NB+i];
            }

    }
}


        double MPIt1 = MPI_Wtime();
        this->LUfactorisation(A_loc,N,ipiv,BWU,BWL,desca,AF,LAF,WORK,LWORK,mycol);
        this->LinearSolver(A_loc,B_loc,N,ipiv,BWU,BWL,desca,descb,AF,LAF,WORK,LWORK,mycol);
        
        double MPIt2 = MPI_Wtime();
        
        double MPIelapsed=MPIt2-MPIt1;
    
        Cblacs_gridexit( ctxt );//done using grid
    
        int rank, psize;
        //---------------------- Look this part again
        MPI_Comm_rank(MPI_COMM_WORLD,&rank);
        MPI_Comm_size(MPI_COMM_WORLD,&psize); 
    
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Gather(B_loc,NB,MPI_DOUBLE,B_loc_red,NB,MPI_DOUBLE,0,MPI_COMM_WORLD);
    
        if (rank==0) {
        cblas_dcopy((Nx-2)*(Ny-2),B_loc_red,1,s_inner,1);
        }
    
        MPI_Barrier(MPI_COMM_WORLD);
        
        Cblacs_exit(ctxt);
    
        if (rank==0) {
        cout<<"Updated s_inner: "<<endl;
        print_matrix(s_inner,Ny-2,Nx-2);
        cout<<"Time taken for PDBTRF/PDGBTRS: "<<endl;
        printf( "%6.20lf \n", MPIelapsed );
        
    }
    
}


void PoissonSolver:: LUfactorisation(double* A_loc, int N, int* ipiv,int BWU,int BWL
        ,int* desca,double* AF, int LAF,double* WORK,int LWORK,int mycol)
{
    int info;
    
    F77NAME(pdgbtrf) (N,BWL,BWU,A_loc,1,desca,ipiv,AF,LAF,WORK,LWORK,&info);
    if (info!=0) {
        cout << "Process: "<<mycol<<"Failed to LU factorise matrix, Info: " <<info<< endl;
    }
}


void PoissonSolver::LinearSolver (double* A_loc,double* B_loc, int N, int* ipiv,int BWU,int BWL
        ,int* desca,int* descb,double* AF, int LAF,double* WORK,int LWORK,int mycol)
{

    int info;
    F77NAME(pdgbtrs) ('N',N,BWL,BWU,1,A_loc,1,
        desca,ipiv,B_loc,1,descb,AF,LAF,WORK,LWORK,&info);
    if (info!=0) {
        cout <<"Process: "<<mycol<< "Failed to solve, Info: " <<info<< endl;
    }
}


