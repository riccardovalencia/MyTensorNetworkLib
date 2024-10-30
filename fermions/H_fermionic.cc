#include "H_fermionic.h"
#include <itensor/all.h>
#include <iostream>


MPO
H_tight_binding_electrons(const int N , const SiteSet sites, const vector<double> J, const vector<double> hup, const vector<double> hdn)
{
    int size_J = J.size();
    int size_hup = hup.size();
    int size_hdn = hdn.size();
    
    
    auto ampo = AutoMPO(sites);

    for(int j : range1(N-1))
    {
        
        ampo += J[(j-1)%size_J] , "Cdagup" , j   , "Cup" , j+1;
        ampo += J[(j-1)%size_J] , "Cdagup" , j+1 , "Cup" , j  ;
        ampo += J[(j-1)%size_J] , "Cdagdn" , j   , "Cdn" , j+1;
        ampo += J[(j-1)%size_J] , "Cdagdn" , j+1 , "Cdn" , j  ;
        
        ampo += hup[(j-1)%size_hup] , "Nup" , j;
        ampo += hdn[(j-1)%size_hdn]  , "Ndn" , j;
        
    }

    // else
    // {

    //     for(int j : range1(N-1))
    //     {
    //         for(int j : range1(N-1))
    //     {
    //         // if(j==1)
    //         // {
    //         //     ampo += V , "Cdag"   , j   , "C", j+1;
    //         //     ampo += V , "Cdag"   , j+1 , "C", j;
                
    //         // }
    //         // else if(j<N)
    //         // {
    //         ampo +=  J , "Cdag" , j   , "C" , j+1;
    //         ampo +=  J , "Cdag" , j+1 , "C" , j;
    //         // }

    //         ampo += 
    //     }    

            
    //     }
    // }


    return toMPO(ampo);
}



