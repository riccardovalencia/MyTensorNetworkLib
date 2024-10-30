#include "initial_state_fermions.h"
#include <itensor/all.h>
#include <iostream>

MPS 
fermi_sea_electrons(MPO H, const SiteSet sites, const int Nupfill, const int Ndnfill, Sweeps sweeps, double min_varH)
{
    int N = length(sites);

    if(Nupfill + Ndnfill > 2*N || Nupfill > N || Ndnfill > N)
    {
        cerr << "Filling larger than the one it can be hosted.\n";
        exit(-1);
    }

    InitState state = InitState(sites,"0");

    for(int j = 1; j <= Nupfill + Ndnfill; j += 1)
    {
        if(j <= Nupfill && j <= Ndnfill) state.set(j,"UpDn");
        else if(j<=Nupfill) state.set(j,"Up");
        else if(j<=Ndnfill) state.set(j,"Dn");
    }

    MPS psi0 = MPS(state); 

    cerr << "Computing ground state...\n";

    auto [energy,psi_gs] = dmrg(H,psi0,sweeps,{"Quiet",true});
    double E2 = inner(H,psi_gs,H,psi_gs) ;
    double E  = inner(psi_gs,H,psi_gs);
    double varE = E2 - E*E;
    
    if(abs(varE) > min_varH)
    {
        cerr << "Variance too large : " << varE << "\n";
        exit(-1);
    }

    // check filling
    AutoMPO ampo = AutoMPO(sites);
    for(int j : range1(N)) ampo += 1, "Nup" , j;
    MPO Nuptot = toMPO(ampo);

    ampo = AutoMPO(sites);
    for(int j : range1(N)) ampo += 1, "Ndn" , j;
    MPO Ndntot = toMPO(ampo);

    double nuptot = real(innerC(psi_gs,Nuptot,psi_gs));
    double ndntot = real(innerC(psi_gs,Ndntot,psi_gs));

     if( abs(Nupfill - nuptot) > 1E-7 ||  abs(Ndnfill - ndntot) > 1E-7 )
    {
        cerr << "Expected : " << Nupfill << " Computed : " << nuptot << "\n";
        cerr << "Expected : " << Ndnfill << " Computed : " << ndntot << "\n";
        exit(-1);
    }


    return psi_gs;

}


MPS
initial_computational_electron_state(const SiteSet sites, const int Nupfill, const int Ndnfill)
{
    int N = length(sites);

    if(Nupfill + Ndnfill > 2*N || Nupfill > N || Ndnfill > N)
    {
        cerr << "Filling larger than the one it can be hosted.\n";
        exit(-1);
    }

    InitState state = InitState(sites,"0");

    for(int j = 1; j <= Nupfill + Ndnfill; j += 1)
    {
        if(j <= Nupfill && j <= Ndnfill) state.set(j,"UpDn");
        else if(j<=Nupfill) state.set(j,"Up");
        else if(j<=Ndnfill) state.set(j,"Dn");
    }

    return MPS(state); 

}
