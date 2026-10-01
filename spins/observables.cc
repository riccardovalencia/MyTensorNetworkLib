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

//measure of longitudinal and trasnversal magnetization in each site

void 
measure_mx_mz( const SpinHalf sites , MPS psi , const int N)
	{
	
	for( int j = 1 ; j <= N ; j++ )
		{
		psi.position(j);
		double Mx1 = 2 * eltC(dag(prime(psi(j),"Site")) * op(sites,"Sx",j) * psi(j)).real();
		double Mz1 = 2 * eltC(dag(prime(psi(j),"Site")) * op(sites,"Sz",j) * psi(j)).real();
		cout << "Sx_" << j << " = " << Mx1 << "\n"
			 << "Sz_" << j << " = " << Mz1 << endl;
		}
	}
