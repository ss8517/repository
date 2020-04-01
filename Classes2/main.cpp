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

using namespace std;

#include "LidDrivenCavity.h"

int main() {   
    double Lx=1.0;
    double Ly=1.0;
    double Re=100;
    int Nx=12;//columns
    int Ny=12;//rows
    double dt=.0001;//time step 
    double T=100;
    int Px=3; 
    int Py=3; 
    
    LidDrivenCavity* solver=new LidDrivenCavity((Nx-2)*(Ny-2),Nx,Ny); 

    solver->Initialise(Lx,Ly,Nx,Ny,Px,Py,dt,T,Re);

    //delete [] solver;
}