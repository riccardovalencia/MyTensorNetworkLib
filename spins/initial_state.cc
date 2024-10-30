#include "initial_state.h"
#include <itensor/all.h>
#include <iostream>

using namespace std;
using namespace itensor;

//----------------------------------------------------------------------
// return MPS representing all spins up along x

void 
initial_state_all_UP( const SpinHalf sites , MPS* psi , const int N)
	{
	for(int i=1; i<=N; i++)
		{
		auto si = sites(i);
		auto wf = ITensor(si);
		
		wf.set(si(1), 1/sqrt(2));
		wf.set(si(2), 1/sqrt(2));
	
		(*psi).setA(i,wf);
		}
	}
	
//----------------------------------------------------------------------
// return MPS representing all spins down along x

void 
initial_state_all_DOWN( const SpinHalf sites , MPS* psi , const int N)
	{
	for(int i=1; i<=N; i++)
		{
		auto si = sites(i);
		auto wf = ITensor(si);
		
		wf.set(si(1), 1/sqrt(2));
		wf.set(si(2), -1/sqrt(2));
	
		(*psi).setA(i,wf);
		}
	}
	
//----------------------------------------------------------------------
// return MPS with half chain DOWN along X, half chain UP along X (domain wall with single kink)

void
initial_state_DOMAIN_WALL( const SpinHalf sites , MPS* psi , const int N)
	{
	for(int i=1; i<=N; i++)
		{
		auto si = sites(i);
		auto wf = ITensor(si);
		if( i <= N/2 )
			{
			wf.set(si(1), 1/sqrt(2));
			wf.set(si(2), 1/sqrt(2));
			}
		else
			{
			wf.set(si(1), 1/sqrt(2));
			wf.set(si(2), -1/sqrt(2));
			}
		(*psi).setA(i,wf);
		}
	}
	


// MPS
// initial_computational_state(const SiteSet sites , const vector<int> initial_state)
// {

//     MPS psi = randomMPS(sites);
//     int N   = length(sites);

//     if(N != initial_state.size())
//     {
//         cerr << "Lenght of initial state is different from N!\nReturning random MPS.\n";
//         return psi;
//     }

//     // for(int q : initial_state) cerr << q << " ";
//     // exit(0);


//     if(N==1)
//     {
//         Index sj = sites(1);
//         ITensor wf = ITensor(sj);
        
        
        
//         if(hasTags(sj,"Site,S=1/2"))
//         {
//             cerr << "Inserting spin coherent state" << endl;
//             if (initial_state[0]==0)
//             {
//                 wf.set(sj(1),1);
//                 wf.set(sj(2),0);
//             }
//             else
//             {
//                 wf.set(sj(1),0);
//                 wf.set(sj(2),1);
//             }
            
//         }
        
//         else{
//             cerr << "SiteSet not recognize : " << sj << endl;
//             cerr << "Return a random initial state" << endl;
//             return psi; 
//         }

//         psi.set(1,wf);
//     }

//     if(N>1)
//     {

//         Index sj ,rj, lj;

//         // first site
//         sj = sites(1);
//         rj = commonIndex(psi(1),psi(2));
//         ITensor wf = ITensor(sj,rj);
        
//         if(hasTags(sj,"Site,S=1/2"))
//         {
//             if (initial_state[0]==0)
//             {
//                 wf.set(sj(1),rj(1),1);
//                 wf.set(sj(2),rj(1),0);
//             }
//             else
//             {
//                 wf.set(sj(1),rj(1),0);
//                 wf.set(sj(2),rj(1),1);
//             }

//         }
        
//         else{
//             cerr << "SiteSet not recognize : " << sj << endl;
//             cerr << "Return a random initial state" << endl;
//             return psi; 
//         }

//         psi.set(1,wf);
        
//         for(int j=2 ; j < N ; j++)
//         {
//             sj = sites(j);
//             lj = commonIndex(psi(j-1),psi(j));
//             rj = commonIndex(psi(j),psi(j+1));
//             wf = ITensor(sj,lj,rj);

//             if(hasTags(sj,"Site,S=1/2"))
//             {
//                 if (initial_state[j-1]==0)
//                 {
//                     cerr << "Inserted spin down\n";
//                     wf.set(sj(1),lj(1),rj(1),1);
//                     wf.set(sj(2),lj(1),rj(1),0);
//                 }
//                 else
//                 {
//                     cerr << "Inserted spin up\n";
//                     wf.set(sj(1),lj(1),rj(1),0);
//                     wf.set(sj(2),lj(1),rj(1),1);
//                 }

//             }
        
//             else
//             {
//                 cerr << "SiteSet not recognize : " << sj << endl;
//                 cerr << "Return a random initial state" << endl;
//                 return psi; 
//             }
                
//             psi.set(j,wf); 
//         }

//         sj = sites(N);
//         lj = commonIndex(psi(N-1),psi(N));
//         wf = ITensor(sj,lj);

                        
//         if(hasTags(sj,"Site,S=1/2"))
//         {
//             if (initial_state[N-1]==0)
//             {
//                 wf.set(sj(1),lj(1),1);
//                 wf.set(sj(2),lj(1),0);
//             }
//             else
//             {
//                 wf.set(sj(1),lj(1),0);
//                 wf.set(sj(2),lj(1),1);
//             }
//         }

//         psi.set(N,wf); 
//     }

//     return psi;
// }



	

	

	
