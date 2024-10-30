#include "build_hamiltonian.h"
#include <itensor/all.h>

using namespace itensor;
using namespace std;



//----------------------------------------------------------------------
// Hamiltonian H = - 0.5 n_0 ( exp(-s) \sigma_1^x -1) - 0.5 \sum_{j=1}^{L-1} n_j (exp(-s) \sigma_{j+1}^x - (1-2c) n_{j+1} -1) + 0.5 n_L * sector + 0.5 * n_L

MPO
H_mmGcbQEM( const SiteSet sites, int size , int n0, double symmetry , double s, double c)
{

	double U = 1-2*c;

	auto ampo = AutoMPO(sites);

	ampo += - n0 * exp(-s) * 0.5  , "A" , 1 ;
	ampo += - n0 * exp(-s) * 0.5  , "Adag" , 1 ;
	ampo += 0.5 * n0 * U , "N", 1 ;
	ampo += n0 * 0.5 , "Id", 1;

	for(int j = 1 ; j <= size-1 ; j++)
		{
		ampo += - exp(-s) * 0.5 , "N" , j , "A" , j+1;
		ampo += - exp(-s) * 0.5 , "N" , j , "Adag" , j+1;
		ampo += 0.5 * U , "N", j , "N" , j+1 ;
		ampo += 0.5 , "N", j , "Id", j+1;
		}

	// eigenvalues for the truly bosonic system

	// double eigenvalue_symmetry;
	// if( U > 0) eigenvalue_symmetry = symmetry * U - exp(-2*s)/U;
	// else eigenvalue_symmetry = 0;

	// eigenvalues for the finite cut-off system
	
	// double eigenvalue_symmetry = symmetry;
	// test
	double eigenvalue_symmetry = -1*symmetry;


	ampo += 0.5 * eigenvalue_symmetry , "N", size; 
	ampo += 0.5, "N", size; 

	MPO H = toMPO(ampo,{"Exact=",true});

	return H;
}


//----------------------------------------------------------------------
// Hamiltonian H = - 0.5 n_0 ( exp(-s) \sigma_1^x -1) - 0.5 \sum_{j=1}^{L-1} n_j (exp(-s) \sigma_{j+1}^x - (1-2c) n_{j+1} -1) + 0.5 n_L * sector + 0.5 * n_L

MPO
H_mmGcbQEM_with_drift( const SiteSet sites, int size , int n0, double symmetry , double s, double c)
{

	double U = 1-2*c;
	double Omega = 0.05;
	auto ampo = AutoMPO(sites);

	ampo += - n0 * exp(-s) * 0.5  , "A" , 1 ;
	ampo += - n0 * exp(-s) * 0.5  , "Adag" , 1 ;
	ampo += 0.5 * n0 * U , "N", 1 ;
	ampo += n0 * 0.5 , "Id", 1;

	for(int j = 1 ; j <= size-1 ; j++)
		{
		ampo += - exp(-s) * 0.5 , "N" , j , "A" , j+1;
		ampo += - exp(-s) * 0.5 , "N" , j , "Adag" , j+1;
		ampo += 0.5 * U , "N", j , "N" , j+1 ;
		ampo += 0.5 , "N", j , "Id", j+1;
		}

	for(int j = 2 ; j <= size-1 ; j++)
		{
		ampo += Omega , "A" , j ;
		ampo += Omega , "Adag" , j;
		}

	// eigenvalues for the truly bosonic system

	// double eigenvalue_symmetry;
	// if( U > 0) eigenvalue_symmetry = symmetry * U - exp(-2*s)/U;
	// else eigenvalue_symmetry = 0;

	// eigenvalues for the finite cut-off system
	
	// double eigenvalue_symmetry = symmetry;
	// test
	double eigenvalue_symmetry = -1*symmetry;


	ampo += 0.5 * eigenvalue_symmetry , "N", size; 
	ampo += 0.5, "N", size; 

	MPO H = toMPO(ampo,{"Exact=",true});

	return H;
}



MPO
H_mmGcbQEM_minus( const SiteSet sites, int size , int n0, double symmetry , double s, double c)
{

	 
	double U = 1-2*c;

	auto ampo_minus = AutoMPO(sites);

	ampo_minus += n0 * exp(-s) * 0.5  , "A" , 1 ;
	ampo_minus += n0 * exp(-s) * 0.5  , "Adag" , 1 ;
	ampo_minus += -0.5 * n0 * U , "N", 1 ;
	ampo_minus += -n0 * 0.5 , "Id", 1;

	for(int j = 1 ; j <= size-1 ; j++)
		{
		ampo_minus += exp(-s) * 0.5 , "N" , j , "A" , j+1;
		ampo_minus += exp(-s) * 0.5 , "N" , j , "Adag" , j+1;
		ampo_minus += -0.5 * U , "N", j , "N" , j+1 ;
		ampo_minus += -0.5 , "N", j , "Id", j+1;
		}

	
	// eigenvalues for the truly bosonic system

	// double eigenvalue_symmetry;
	// if( U > 0) eigenvalue_symmetry = symmetry * U - exp(-2*s)/U;
	// else eigenvalue_symmetry = 0;

	// eigenvalues for the finite cut-off system
	
	double eigenvalue_symmetry = symmetry;

	ampo_minus += -0.5 * eigenvalue_symmetry , "N", size; 
	ampo_minus += -0.5, "N", size;

	MPO H_minus = toMPO(ampo_minus,{"Exact=",true});

	return H_minus;

}


MPO
H_mmGcbQEM_onsite_hopping( const SiteSet sites, int size , double s, double c, double epsilon, double t)
{

	 
	double U = 1-2*c;
	auto ampo = AutoMPO(sites);


	for(int j = 1 ; j <= size-1 ; j++)
		{
		ampo += - exp(-s) * 0.5 , "N" , j , "A" , j+1;
		ampo += - exp(-s) * 0.5 , "N" , j , "Adag" , j+1;
		ampo += 0.5 * U , "N", j , "N" , j+1 ;
		ampo += 0.5 , "N", j , "Id", j+1;

		// on site density-density
		ampo += 0.5 * epsilon , "N", j , "N", j;

		ampo += -t * 0.5, "Adag", j , "A", j+1;
		ampo += -t * 0.5, "A", j+1 , "Adag", j;

		}

	// on site density-density
	ampo += 0.5 * epsilon , "N", size , "N", size;

	MPO H = toMPO(ampo,{"Exact=",true});

	return H;
}


MPO
H_mmGcbQEM_onsite( const SiteSet sites, int size , int n0, double symmetry , double s, double c , double epsilon)
{
	 
	double U = 1-2*c;
	auto ampo = AutoMPO(sites);

	// ampo += - n0 * exp(-s) * 0.5  , "A" , 1 ;
	// ampo += - n0 * exp(-s) * 0.5  , "Adag" , 1 ;
	// ampo += 0.5 * n0 * U , "N", 1 ;
	// ampo += n0 * 0.5 , "Id", 1;

	for(int j = 1 ; j <= size-1 ; j++)
		{
		// ampo += - exp(-s) * 0.5 , "N" , j , "A" , j+1;
		// ampo += - exp(-s) * 0.5 , "N" , j , "Adag" , j+1;
		// ampo += 0.5 * U , "N", j , "N" , j+1 ;
		// ampo += 0.5 , "N", j ;
		ampo += 0.5 * epsilon, "N", j , "N", j ;
		}

	ampo += 0.5 * epsilon , "N", size , "N", size;
	// ampo += 0.5, "N", size;	

	MPO H = toMPO(ampo,{"Exact=",true});

	return H;
}


// MPO
// H_mmGcbQEM_onsite( const SiteSet sites, int size , int n0, double symmetry , double s, double c , double epsilon)
// {
	 
// 	double U = 1-2*c;
// 	auto ampo = AutoMPO(sites);

// 	ampo += - n0 * exp(-s) * 0.5  , "A" , 1 ;
// 	ampo += - n0 * exp(-s) * 0.5  , "Adag" , 1 ;
// 	ampo += 0.5 * n0 * U , "N", 1 ;
// 	ampo += n0 * 0.5 , "Id", 1;

// 	// non serve mettere n_0^2 dato che e' una costante.
// 	// epsilon = 2;

// 	for(int j = 1 ; j <= size-1 ; j++)
// 		{
// 		ampo += - exp(-s) * 0.5 , "N" , j , "A" , j+1;
// 		ampo += - exp(-s) * 0.5 , "N" , j , "Adag" , j+1;
// 		ampo += 0.5 * U , "N", j , "N" , j+1 ;
// 		ampo += 0.5 , "N", j ;
// 		ampo += 0.5 * epsilon, "N", j , "N", j ;
// 		}

// 	// on site density-density
// 	ampo += 0.5 * epsilon , "N", size , "N", size;
// 	// ampo += 0.5 * epsilon , "N", size;

// 	// exit(0);
// 	// la presenza di on-site density-density interaction non cambia il settore della simmetria del modello che non lo ha. La popolazione sul sito L+1 e' solo un offset
// 	// in questo caso, dato che non parla con n_L. 

	
// 	// eigenvalues for the truly bosonic system

// 	// double eigenvalue_symmetry;
// 	// if( U > 0) eigenvalue_symmetry = symmetry * U - exp(-2*s)/U;
// 	// else eigenvalue_symmetry = 0;

// 	// eigenvalues for the finite cut-off system
	
// 	double eigenvalue_symmetry = symmetry;

// 	// ampo += 0.5 * eigenvalue_symmetry , "N", size; 
// 	ampo += 0.5, "N", size;	

// 	MPO H = toMPO(ampo,{"Exact=",true});

// 	return H;
// }


MPO
H_mmGcbQEM_onsite_nonext( const SiteSet sites, int size , int n0, double symmetry , double s, double c )
{

	auto ampo = AutoMPO(sites);

	// ampo += - n0 * exp(-s) * 0.5  , "A" , 1 ;
	// ampo += - n0 * exp(-s) * 0.5  , "Adag" , 1 ;
	ampo +=   n0 * 0.5 , "Id", 1;

	// non serve mettere n_0^2 dato che e' una costante.

	for(int j = 1 ; j <= size-1 ; j++)
		{
		// ampo += - exp(-s) * 0.5 , "N" , j , "A" , j+1;
		// ampo += - exp(-s) * 0.5 , "N" , j , "Adag" , j+1;
		ampo += 0.5 , "N", j , "Id", j+1;

		// on site density-density
		ampo += 0.5 * c , "N", j , "N", j;

		}

	// on site density-density
	ampo += 0.5 * c , "N", size , "N", size;
	ampo += - symmetry * 0.5 , "N", size; 
	ampo += 0.5, "N", size;

	MPO H = toMPO(ampo,{"Exact=",true});


	return H;
}

// ------------------------------------------------------------
// Build Hamiltonian of Bosonic quantum east model with next-neighbour density-density interaction.
// No on-site density-density interaction
// site 0 is not touched. This is used to find the ground state within a certain symmetry sector fixed by
// n0 from a site of the form |n0>|psi>.

MPO
H_mmGcbQEM_n0_untouched( const SiteSet sites, int size , int n0, double symmetry , double s, double c)
{

	 
	double U  = 1 - 2 * c;
	auto ampo = AutoMPO(sites);

	ampo += - n0 * exp(-s) * 0.5  , "A" , 2 ;
	ampo += - n0 * exp(-s) * 0.5  , "Adag" , 2 ;
	ampo +=  0.5 * n0 * U , "N", 2 ;
	ampo +=   n0 * 0.5 , "Id", 2;

	for(int j = 2 ; j <= size-1 ; j++)
		{
		ampo += - 0.5 * exp(-s) , "N" , j , "A"    , j+1;
		ampo += - 0.5 * exp(-s) , "N" , j , "Adag" , j+1;
		ampo +=   0.5 * U , "N", j , "N" , j+1 ;
		ampo +=   0.5     , "N", j ;
		}

	
	// eigenvalues for the truly bosonic system

	// double eigenvalue_symmetry;
	// if( U > 0) eigenvalue_symmetry = symmetry * U - exp(-2*s)/U;
	// else eigenvalue_symmetry = 0;

	// eigenvalues for the finite cut-off system
	
	double eigenvalue_symmetry = -1*symmetry;
 
	ampo += 0.5 * eigenvalue_symmetry , "N", size; 
	ampo += 0.5, "N", size;

	// ampo += - symmetry , "N", size; 
	// ampo +=   0.5     , "N", size ;

	MPO H = toMPO(ampo,{"Exact=",true});

	return H;
}

// ------------------------------------------------------------
// Build Hamiltonian of Bosonic quantum east model with next-neighbour density-density interaction.
// No on-site density-density interaction
// Notice that n0 is not a number, but its defined in the Hamiltonian itself. So, it is not a conserved quantity
// and we do not expect DMRG to conserve it (it will try to minimize the Hamiltonian and that's it. At most it 
// will be stuck for a long time if you initialize DMRG with a state |n_0> \otimes |\psi> since |n_0> is an eigenstate
// of the operators acting on the 0-th site). 

MPO
H_mmGcbQEM_n0_not_fixed( const SiteSet sites, int size , double symmetry , double s, double c)
{

	 
	double U = 1-2*c;
	auto ampo = AutoMPO(sites);

	for(int j = 1 ; j <= size-1 ; j++)
		{
		ampo += - 0.5 * exp(-s) , "N" , j , "A"    , j+1;
		ampo += - 0.5 * exp(-s) , "N" , j , "Adag" , j+1;
		ampo +=   0.5 * U , "N", j , "N" , j+1 ;
		ampo +=   0.5     , "N", j ;
		}


	
	// eigenvalues for the truly bosonic system

	// double eigenvalue_symmetry;
	// if( U > 0) eigenvalue_symmetry = symmetry * U - exp(-2*s)/U;
	// else eigenvalue_symmetry = 0;

	// eigenvalues for the finite cut-off system
	
	double eigenvalue_symmetry = -1*symmetry;

	ampo += 0.5 * eigenvalue_symmetry , "N", size; 
	ampo += 0.5, "N", size;

	// ampo += - symmetry , "N", size; 
	// ampo +=   0.5     , "N", size ;

	MPO H = toMPO(ampo,{"Exact=",true});

	return H;
}


// Exponential of the Hamiltonian. This is done in order to apply it to operators

MPO
exp_H_mmGcbQEM_n0_notfixed( const SiteSet sites, int size , int n0, double symmetry , double J, double c, double dt)
{
	 
	double U  = 1 - 2 * c;
	auto ampo = AutoMPO(sites);

	for(int j = 1 ; j <= size-1 ; j++)
		{
		ampo += - 0.5 * J , "N" , j , "A"    , j+1;
		ampo += - 0.5 * J , "N" , j , "Adag" , j+1;
		ampo +=   0.5 * U , "N", j , "N" , j+1 ;
		ampo +=   0.5     , "N", j ;
		}

	// eigenvalues for the truly bosonic system

	// double eigenvalue_symmetry;
	// if( U > 0) eigenvalue_symmetry = symmetry * U - exp(-2*s)/U;
	// else eigenvalue_symmetry = 0;

	// eigenvalues for the finite cut-off system
	
	double eigenvalue_symmetry = -1*symmetry;

	ampo += 0.5 * eigenvalue_symmetry , "N", size; 
	ampo += 0.5, "N", size;	

	MPO expH = toExpH( ampo , Cplx_i * dt );
	return expH;
}