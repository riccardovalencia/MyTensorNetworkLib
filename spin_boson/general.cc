#include "general.h"
#include <itensor/all.h>
using namespace std;
using namespace itensor;


// ----------------------------------------------------------
// Apply a gate on the jn sites of an MPS psi 
// It is possible to apply up to 3-sites gates.

MPS
apply_gate(MPS psi, const ITensor gate, const vector<int> jn, const Args args)
{
    double cut_off = args.getReal("Cutoff");
    int maxDim  = args.getInt("MaxDim");
	int j = jn[0];

	if( j < 0 || j > length(psi))
	{
		cerr << "Site not valid. Cannot apply gate \n";
		exit(0);
	}

	ITensor AA = gate;
	psi.position(j);
           
	for(int q : jn) AA *= psi(q);  
	AA.mapPrime(1,0);

	if(jn.size() == 1)
	{
		psi.set(j,AA);
	}

	else if(jn.size() == 2)
	{
		auto [U,S,V] = svd(AA,inds(psi(j)),{"Cutoff=",cut_off,"MaxDim=",maxDim});
		psi.set(j,U);
		psi.set(j+1,S*V);
	}

	else if(jn.size() == 3)
	{
		auto [U,S,V] = svd(AA,inds(psi(j)),{"Cutoff=",cut_off,"MaxDim=",maxDim});
		Index l =  commonIndex(U,S);
		psi.set(j,U);

		auto [U1,S1,V1] = svd(S*V,{inds(psi(j+1)),l},{"Cutoff=",cut_off,"MaxDim=",maxDim});
		psi.set(j+1,U1);
		psi.set(j+2,S1*V1);
	}

	else
	{
	cerr << jn.size() <<"-gate not implemented (yet)!\n";
	exit(0);
	}
    
	return psi;
}