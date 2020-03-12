#include "LidDrivenCavity.h"

LidDrivenCavity::LidDrivenCavity()
{
}

LidDrivenCavity::~LidDrivenCavity()
{
}

void LidDrivenCavity::SetDomainSize(double xlen, double ylen)
{
    Lx=xlen;
    Ly=ylen;
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

void LidDrivenCavity::Initialise()
{
    
}

void LidDrivenCavity::Integrate()
{
}
