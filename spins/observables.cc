#include "observables.h"
#include <itensor/all.h>

#include <sys/stat.h>
#include <iostream>
#include <string>
#include <vector>
#include <ctime>
#include <fstream>	//output file
#include <sstream>	//for ostringstream
#include <iomanip>

using namespace itensor;
using namespace std;

//----------------------------------------------------------------------

//measure of entanglement entropy centered in a certain site

double
entanglement_entropy( MPS* psi , int N , int site)
	{
	(*psi).position(site); 
	ITensor wf = (*psi).A(site) * (*psi).A(site+1);
	ITensor U  = (*psi).A(site);
	ITensor S,V;
	auto spectrum = svd(wf,U,S,V);
	
	double SvN = 0.;
	for(auto p : spectrum.eigs())
		{
		if(p > 1E-12) SvN += -p*log(p);
		}
	return SvN;
	}

//----------------------------------------------------------------------

//measure of longitudinal and trasnversal magnetization in each site

void 
measure_mx_mz( const SpinHalf sites , MPS psi , const int N)
	{
	
	for( int j = 1 ; j <= N ; j++ )
		{
		psi.position(j);
		Real Mx1 = 2 * (dag(prime(psi.A(j), Site )) * sites.op("Sx",j) * psi.A(j)).real();
		Real Mz1 = 2 * (dag(prime(psi.A(j), Site )) * sites.op("Sz",j) * psi.A(j)).real();
		cout << "Sx_" << j << " = " << Mx1 << "\n"
			 << "Sz_" << j << " = " << Mz1 << endl;
		}
	}
