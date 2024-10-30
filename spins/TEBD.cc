#include "TEBD.h"
#include <itensor/all.h>
#include <iostream>

using namespace std;
using namespace itensor;

//----------------------------------------------------------------------

//single gate acting on sites [b,b+1] of the Ising model with longitudinal (hx) and transversal (hz) magnetic fields

void
build_single_step( ITensor *hterm , const SpinHalf sites , const int N , const double J , const double hx , const double hz , const int b )
	{
	ITensor Sx1 = sites.op("Sx",b);					
	ITensor Sx2 = sites.op("Sx",b+1);
	ITensor Sz1 = sites.op("Sz",b);
	ITensor Sz2 = sites.op("Sz",b+1);
	ITensor Id1 = sites.op("Id",b);
	ITensor Id2 = sites.op("Id",b+1);
		
	*hterm = - 4 * J * Sx1 * Sx2;
		
	if( b == 1 )
		{
		*hterm +=  - 2 * J * hx * ( Sx1 * Id2 + Id1 * Sx2 / 2. ); 									
		*hterm +=  - 2 * J * hz * ( Sz1 * Id2 + Id1 * Sz2 / 2. );	
		}
	else if( b == N-1)
		{
		*hterm +=  - 2 * J * hx * ( Sx1 * Id2 / 2. + Id1 * Sx2 ); 									
		*hterm +=  - 2 * J * hz * ( Sz1 * Id2 / 2. + Id1 * Sz2 );		
		}
	else{
		*hterm +=  - 2 * J * hx * ( Sx1 * Id2 + Id1 * Sx2 ) / 2.; 									
		*hterm +=  - 2 * J * hz * ( Sz1 * Id2 + Id1 * Sz2 ) / 2.;	
		}
	}
