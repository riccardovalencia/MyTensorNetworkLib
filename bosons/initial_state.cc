#include "initial_state.h"
#include "TEBD.h"
#include "build_hamiltonian.h"
#include "perform_dmrg.h"
#include <itensor/all.h>
#include <iostream>
#include <math.h>

using namespace std;
using namespace itensor;

//----------------------------------------------------------------------
//all spins UP along X

void 
initial_state_all_UP( const SpinHalf sites , MPS* psi , const int N)
	{
	for(int i=1; i<=N; i++)
		{
		auto si = sites(i);
		auto wf = ITensor(si);
		
		wf.set(si(1), 1/sqrt(2));
		wf.set(si(2), 1/sqrt(2));
	
		(*psi).set(i,wf);
		}
	}
	
//----------------------------------------------------------------------
//all spins DOWN along X

void 
initial_state_all_DOWN( const SpinHalf sites , MPS* psi , const int N)
	{
	for(int i=1; i<=N; i++)
		{
		auto si = sites(i);
		auto wf = ITensor(si);
		
		wf.set(si(1), 1/sqrt(2));
		wf.set(si(2), -1/sqrt(2));
	
		(*psi).set(i,wf);
		}
	}
	
//----------------------------------------------------------------------
//half chain DOWN along X, half chain UP along X (domain wall with single kink)
	
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
		(*psi).set(i,wf);
		}
	}


// ----------------------------------------------------------------------
// Initial bosonic state of the form |n0> |0000...0>.
void
initial_state_n0_excitation( MPS* psi, const SiteSet sites, const int size , const int n0 , const int excitation_position)
{
	
	// SITE 1

	Index sj = sites(1);
	Index li = commonIndex((*psi)(1),(*psi)(2));
	ITensor wf = ITensor(sj,li);

	if(excitation_position == 1)
	{
		for( int d=1; d <= n0; d++) wf.set(sj(d),li(1), 0);
		wf.set(sj(n0+1),li(1),1);
		for( int d=n0+2; d <= dim(sj); d++) wf.set(sj(d),li(1), 0);
	}
	else
	{
		wf.set(sj(1),li(1), 1);
		for( int d=2; d <= dim(sj); d++) wf.set(sj(d),li(1), 0);
	}
	(*psi).set(1,wf);

	// SITE 2 TO SIZE-1

	for( int j = 2 ; j < size; j++)
	{
		sj = sites(j);
		li = commonIndex((*psi)(j-1),(*psi)(j));
		Index ri = commonIndex((*psi)(j),(*psi)(j+1));
		wf = ITensor(sj,li,ri);

		if(j==excitation_position)
		{
			for( int d=1; d <= n0; d++) wf.set(sj(d),li(1),ri(1), 0);
			wf.set(sj(n0+1),li(1),ri(1),1);
			for( int d=n0+2; d <= dim(sj); d++) wf.set(sj(d),li(1),ri(1), 0);
		}
		else
		{
			wf.set(sj(1),li(1),ri(1), 1);
			for( int d=2; d <= dim(sj); d++) wf.set(sj(d),li(1),ri(1), 0);
		}
		(*psi).set(j,wf);
	}

	// SITE SIZE (LAST ONE)
	sj  = sites(size);
	li = commonIndex((*psi)(size-1),(*psi)(size));
	wf = ITensor(sj,li);

	if(excitation_position == size)
	{
		for( int d=1; d <= n0; d++) wf.set(sj(d),li(1), 0);
		wf.set(sj(n0+1),li(1),1);
		for( int d=n0+2; d <= dim(sj); d++) wf.set(sj(d),li(1), 0);
	}
	else
	{
		wf.set(sj(1),li(1), 1);
		for( int d=2; d <= dim(sj); d++) wf.set(sj(d),li(1), 0);
	}
	(*psi).set(size,wf);

}


void
initial_state_n0_excitation_pinned( MPS* psi, const SiteSet sites, const int size , const int n0 , const int excitation_position)
{
	
	// SITE 1

	Index sj = sites(1);
	Index li = commonIndex((*psi)(1),(*psi)(2));
	ITensor wf = ITensor(sj,li);

	if(excitation_position == 1)
	{
		sj = sites(1);
		li = commonIndex((*psi)(1),(*psi)(2));
		wf = ITensor(sj,li);
		for( int d=1; d <= n0; d++) wf.set(sj(d),li(1), 0);
		wf.set(sj(n0+1),li(1),1);
		for( int d=n0+2; d <= dim(sj); d++) wf.set(sj(d),li(1), 0);
		(*psi).set(1,wf);
	}
	

	// SITE 2 TO SIZE-1

	else if(excitation_position > 1 && excitation_position < size)
	{
		sj = sites(excitation_position);
		li = commonIndex((*psi)(excitation_position-1),(*psi)(excitation_position));
		Index ri = commonIndex((*psi)(excitation_position),(*psi)(excitation_position+1));
		wf = ITensor(sj,li,ri);

		
		for( int d=1; d <= n0; d++) wf.set(sj(d),li(1),ri(1), 0);
		wf.set(sj(n0+1),li(1),ri(1),1);
		for( int d=n0+2; d <= dim(sj); d++) wf.set(sj(d),li(1),ri(1), 0);
		
		(*psi).set(excitation_position,wf);
	}

	// SITE SIZE (LAST ONE)
	
	if(excitation_position == size)
	{
		sj  = sites(size);
		li = commonIndex((*psi)(size-1),(*psi)(size));
		wf = ITensor(sj,li);
		for( int d=1; d <= n0; d++) wf.set(sj(d),li(1), 0);
		wf.set(sj(n0+1),li(1),1);
		for( int d=n0+2; d <= dim(sj); d++) wf.set(sj(d),li(1), 0);
		(*psi).set(size,wf);
	}


}


// ----------------------------------------------------------------------
// Initial vacuum bosonic state |0000...0>.
void
initial_state_vacuum_state( MPS* psi, const SiteSet sites, const int size )
{
	for( int j=1 ; j<= size; j++)
	{
		Index sj   = sites(j);
		ITensor wf = ITensor(sj);
		wf.set(sj(1), 1);
		for( int d=2; d <= dim(sj); d++) wf.set(sj(d), 0);
		(*psi).set(j,wf);
	}
}

// ----------------------------------------------------------------------
// Initial vacuum bosonic state |0000...0>.
void
initial_state_vacuum_state_correct_link( MPS* psi, const SiteSet sites, const int size )
{

	Index sj = sites(1);
	Index li = commonIndex((*psi)(1),(*psi)(2));
	ITensor wf = ITensor(sj,li);
	wf.set(sj(1),li(1), 1);
	for( int d=2; d <= dim(sj); d++) wf.set(sj(d),li(1), 0);
	(*psi).set(1,wf);

	for( int j = 2 ; j < size; j++)
	{
		sj = sites(j);
		li = commonIndex((*psi)(j-1),(*psi)(j));
		Index ri = commonIndex((*psi)(j),(*psi)(j+1));
		wf = ITensor(sj,li,ri);
		wf.set(sj(1),li(1),ri(1), 1);
		for( int d=2; d <= dim(sj); d++) wf.set(sj(d),li(1),ri(1), 0);
		(*psi).set(j,wf);
	}

	sj  = sites(size);
	li = commonIndex((*psi)(size-1),(*psi)(size));
	wf = ITensor(sj,li);
	wf.set(sj(1),li(1), 1);
	for( int d=2; d <= dim(sj); d++) wf.set(sj(d),li(1), 0);
	(*psi).set(size,wf);
}


// ----------------------------------------------------------------------
// Initial all one bosonic state |1111...1>.
void
initial_state_all_one_state_correct_link( MPS* psi, const SiteSet sites, const int size )
{

	Index sj = sites(1);
	Index li = commonIndex((*psi)(1),(*psi)(2));
	ITensor wf = ITensor(sj,li);
	wf.set(sj(1),li(1), 0);
	wf.set(sj(2),li(1), 1);
	for( int d=3; d <= dim(sj); d++) wf.set(sj(d),li(1), 0);
	(*psi).set(1,wf);

	for( int j = 2 ; j < size; j++)
	{
		sj = sites(j);
		li = commonIndex((*psi)(j-1),(*psi)(j));
		Index ri = commonIndex((*psi)(j),(*psi)(j+1));
		wf = ITensor(sj,li,ri);
		wf.set(sj(1),li(1),ri(1), 0);
		wf.set(sj(2),li(1),ri(1), 1);
		for( int d=3; d <= dim(sj); d++) wf.set(sj(d),li(1),ri(1), 0);
		(*psi).set(j,wf);
	}

	sj  = sites(size);
	li = commonIndex((*psi)(size-1),(*psi)(size));
	wf = ITensor(sj,li);
	wf.set(sj(1),li(1), 0);
	wf.set(sj(2),li(1), 1);
	for( int d=3; d <= dim(sj); d++) wf.set(sj(d),li(1), 0);
	(*psi).set(size,wf);
}


// Initial bosonic state of the form |000..0> |alpha>_j |0000...0>, such that a|alpha> = alpha|alpha>.
void
coherent_state_site_j( MPS* psi, const SiteSet sites, const int size , const int site, const complex<double> alpha )
{

	Index sj = sites(1);
	Index li = commonIndex((*psi)(1),(*psi)(2));
	ITensor wf = ITensor(sj,li);

	if(site==1)
		{
		for( int d=1; d <= dim(sj); d++) wf.set(sj(d),li(1), weigth_coherent_state(alpha, d-1));
		}
	else
		{
		wf.set(sj(1), li(1), 1);
		for( int d=2; d <= dim(sj); d++) wf.set(sj(d),li(1), 0);
		}
	(*psi).set(1,wf);

	for( int j = 2 ; j < size; j++)
	{
		sj = sites(j);
		li = commonIndex((*psi)(j-1),(*psi)(j));
		Index ri = commonIndex((*psi)(j),(*psi)(j+1));
		wf = ITensor(sj,li,ri);
		if(j==site)
			{
			for( int d=1; d <= dim(sj); d++) wf.set(sj(d), li(1), ri(1), weigth_coherent_state(alpha, d-1));
			}
		else
			{
			wf.set(sj(1),li(1),ri(1), 1);
			for( int d=2; d <= dim(sj); d++) wf.set(sj(d), li(1), ri(1), 0);
			}
		(*psi).set(j,wf);
	}

	sj  = sites(size);
	li = commonIndex((*psi)(size-1),(*psi)(size));
	wf = ITensor(sj,li);
	
	if(site==size)
		{
		for( int d=1; d <= dim(sj); d++) wf.set(sj(d),li(1), weigth_coherent_state(alpha, d-1));
		}
	else
		{
		wf.set(sj(1),li(1), 1);
		for( int d=2; d <= dim(sj); d++) wf.set(sj(d),li(1), 0);
		}
	
	(*psi).set(size,wf);


	// -----------------------
	// for( int j=1 ; j<= size; j++)
	// {
	// 	Index sj   = sites(j);
	// 	ITensor wf = ITensor(sj);
	// 	if(j==site)
	// 	{
	// 		// wf.set(sj(n0+1),1);
	// 		// for( int d=n0+2; d <= dim(sj); d++) wf.set(sj(d), 0);
	// 	}
	// 	else
	// 	{
	// 		wf.set(sj(1), 1);
	// 		for( int d=2; d <= dim(sj); d++) wf.set(sj(d), 0);
	// 	}
	// 	(*psi).set(j,wf);
	// }
}


void
coherent_state_all_sites( MPS* psi, const SiteSet sites, const int size , complex<double> alpha )
{

	Index sj = sites(1);
	Index li = commonIndex((*psi)(1),(*psi)(2));
	ITensor wf = ITensor(sj,li);


	
	for( int d=1; d <= dim(sj); d++) wf.set(sj(d),li(1), weigth_coherent_state(alpha, d-1));
	(*psi).set(1,wf);


	for( int j = 2 ; j < size; j++)
	{
		sj = sites(j);
		li = commonIndex((*psi)(j-1),(*psi)(j));
		Index ri = commonIndex((*psi)(j),(*psi)(j+1));
		wf = ITensor(sj,li,ri);
		// for( int d=1; d <= dim(sj); d++) wf.set(sj(d), li(1), ri(1), weigth_coherent_state(alpha/j, d-1));
		for( int d=1; d <= dim(sj); d++) wf.set(sj(d), li(1), ri(1), weigth_coherent_state(alpha, d-1));
		(*psi).set(j,wf);

	}

	sj  = sites(size);
	li = commonIndex((*psi)(size-1),(*psi)(size));
	wf = ITensor(sj,li);
	
	
	for( int d=1; d <= dim(sj); d++) wf.set(sj(d),li(1), weigth_coherent_state(alpha, d-1));
	(*psi).set(size,wf);
}


// Initial bosonic state of the form |000..0> |alpha>_j |0000...0>, such that a|alpha> = alpha|alpha>.
void
squeezed_state_site_j( MPS* psi, const SiteSet sites, const int size , const int site, const double r )
{

	Index sj = sites(1);
	


	if( size > 1)
	{
		Index li = commonIndex((*psi)(1),(*psi)(2));
		ITensor wf = ITensor(sj,li);
		if(site==1)
			{
			for( int d=1; d <= dim(sj); d+=2) wf.set(sj(d),li(1),  weigth_squeezed_state(r, d-1));
			for( int d=2; d <= dim(sj); d+=2) wf.set(sj(d),li(1),  0);
			}
		else
			{
			wf.set(sj(1),li(1), 1);
			for( int d=2; d <= dim(sj); d++) wf.set(sj(d),li(1), 0);
			}
		(*psi).set(1,wf);

		for( int j = 2 ; j < size; j++)
		{
			sj = sites(j);
			li = commonIndex((*psi)(j-1),(*psi)(j));
			Index ri = commonIndex((*psi)(j),(*psi)(j+1));
			wf = ITensor(sj,li,ri);
			if(j==site)
				{
				for( int d=1; d <= dim(sj); d+=2) wf.set(sj(d),li(1),ri(1), weigth_squeezed_state(r, d-1));
				for( int d=2; d <= dim(sj); d+=2) wf.set(sj(d),li(1),ri(1),  0);
				}
			else
				{
				wf.set(sj(1),li(1),ri(1), 1);
				for( int d=2; d <= dim(sj); d++) wf.set(sj(d),li(1),ri(1), 0);
				}
			(*psi).set(j,wf);
		}

		sj  = sites(size);
		li = commonIndex((*psi)(size-1),(*psi)(size));
		wf = ITensor(sj,li);
		
		if(site==size)
			{
			for( int d=1; d <= dim(sj); d+=2) wf.set(sj(d),li(1), weigth_squeezed_state(r, d-1));
			for( int d=2; d <= dim(sj); d+=2) wf.set(sj(d),li(1),  0);
			}
		else
			{
			wf.set(sj(1),li(1), 1);
			for( int d=2; d <= dim(sj); d++) wf.set(sj(d),li(1), 0);
			}
		
		(*psi).set(size,wf);
	}

	else
	{
		ITensor wf = ITensor(sj);
		for( int d=1; d <= dim(sj); d+=2) wf.set(sj(d),  weigth_squeezed_state(r, d-1));
		for( int d=2; d <= dim(sj); d+=2) wf.set(sj(d),  0);
		(*psi).set(size,wf);
	}


}


void
initial_state_cat_state_site_j( MPS* psi, const SiteSet sites, const int size , const int site, const complex<double> alpha )
{
	MPS psi_t_1 = randomMPS(sites);
	MPS psi_t_2 = randomMPS(sites);
	coherent_state_site_j( &psi_t_1, sites, size , site, alpha );
	coherent_state_site_j( &psi_t_2, sites, size , site, -alpha );
	*psi = sqrt(0.5) * sum(psi_t_1,psi_t_2);
	(*psi).position(1);
	(*psi).normalize();
}


// Insert a state within another state, such that you have a state |state_to_insert> that you want to put in another state |psi_t0> from site start to start+L

void
insert_state(MPS* psi_t0, MPS state_to_insert, const SiteSet sites, const SiteSet sites_state_to_insert, const int start, const int L, const int N)
{
    vector<ITensor> copy_of_state_to_insert; 


    Index leftindexj ;
    Index rightindexj ;
    Index physical ;
	Index physical_state_to_insert;
    ITensor psi_tocopy_j ;
    ITensor Tj ;

	int L_effective = 0;

	for(int j=1 ; j<=L ; j++)
	{
		L_effective += 1;
		int current_position_psi_t0 = start + j -1;
		if(current_position_psi_t0<N)
		{
			if(j==1)     leftindexj  = leftLinkIndex(  *psi_t0, current_position_psi_t0);
			else         leftindexj  = leftLinkIndex(  state_to_insert, j);
			if(j==L)     rightindexj = rightLinkIndex( *psi_t0, current_position_psi_t0);
			else         rightindexj = rightLinkIndex( state_to_insert, j);
			physical    = sites(current_position_psi_t0);
			physical_state_to_insert    = sites_state_to_insert(j);
			psi_tocopy_j = state_to_insert(j);

			Tj = ITensor(leftindexj, physical, rightindexj);
			
			for(int l=1; l<= dim(leftindexj) ; l++)
			{
				for(int r=1; r<= dim(rightindexj); r++)
				{
					for(int d=1; d<= dim(physical); d++)
					{
						if(j!=1 && j!=L) Tj.set(leftindexj=l,physical=d,rightindexj=r , eltC(psi_tocopy_j, leftindexj=l,physical_state_to_insert=d, rightindexj=r) );
						else if(j==1) Tj.set(leftindexj=l,physical=d,rightindexj=r , eltC(psi_tocopy_j, physical_state_to_insert=d, rightindexj=r) );			
						else if(j==L) Tj.set(leftindexj=l,physical=d,rightindexj=r , eltC(psi_tocopy_j, physical_state_to_insert=d, leftindexj=l) );
					}
				}
			}

			copy_of_state_to_insert.push_back(Tj);
		}
		
		else if(current_position_psi_t0==N)
		{
			cerr << "Siamo a : " << current_position_psi_t0 << endl;
			leftindexj  = leftLinkIndex( state_to_insert, j);
			rightindexj  = rightLinkIndex( state_to_insert, j);
			physical = sites(current_position_psi_t0);
			physical_state_to_insert = sites_state_to_insert(j);
			psi_tocopy_j = state_to_insert(j);
			cerr << psi_tocopy_j << endl;
			Tj = ITensor(leftindexj, physical);
			
			for(int l=1; l<= dim(leftindexj) ; l++)
			{
					for(int d=1; d<= dim(physical); d++)
					{
						Tj.set(leftindexj=l,physical=d , eltC(psi_tocopy_j, leftindexj=l,physical_state_to_insert=d, rightindexj=1));					
					}
			}

			copy_of_state_to_insert.push_back(Tj);
		}

		else break;
	}
	cerr << "L effective : " << L_effective << endl;	
			
    for(int j=1 ; j<=L_effective ; j++) (*psi_t0).set(j+start-1, copy_of_state_to_insert[j-1]);	

}

// Initial state |0>^k \otimes |n_0> \otimes |GS(n_0)_L> \otimes |0>^(N-L-k-1)
// It's a state with k 0's, n_0, the ground state of size L in this symmetry sector, and the other 0's.
// The total size of the system is L

MPS
super_bosonic_state(MPS state_to_insert, const SiteSet sites, const SiteSet sites_state_to_insert, const int L, const int N, const int n0, const int k)
{
	MPS psi_t0 = randomMPS(sites);
	initial_state_n0_excitation( &psi_t0, sites, N , n0, k+1);
	insert_state(&psi_t0, state_to_insert, sites, sites_state_to_insert, k+2, L, N);
	return psi_t0;
}


// void
// insert_0s_in_the_state(MPS *psi, const int cut_off_fock_space, const int size, const int initial_site, const int final_site)
// {

// 	for(int j = initial_site; j<= final_site; j++)
// 	{
// 		Index sj = sites(j);
// 		if(j==1)
// 		{
// 			Index li = commonIndex((*psi)(1),(*psi)(2));
// 			ITensor wf = ITensor(sj,li);
// 			wf.set(sj(1),li(1), 1);
// 			for( int d=2; d <= cut_off_fock_space; d++) wf.set(sj(d),li(1), 0);
// 			(*psi).set(1,wf);
// 		}
// 		else if(j != size)
// 		{
// 			Index li = commonIndex((*psi)(j-1),(*psi)(j));
// 			Index ri = commonIndex((*psi)(j),(*psi)(j+1));
// 			ITensor wf = ITensor(sj,li,ri);
// 			wf.set(sj(1),li(1),ri(1), 1);
// 			for( int d=2; d <= cut_off_fock_space; d++) wf.set(sj(d),li(1),ri(1), 0);
// 			(*psi).set(j,wf);
// 		}
// 		else if(j==size)
// 		{
// 			Index li = commonIndex((*psi)(size-1),(*psi)(size));
// 			ITensor wf = ITensor(sj,li);
// 			wf.set(sj(1),li(1), 1);
// 			for( int d=2; d <= cut_off_fock_space; d++) wf.set(sj(d),li(1), 0);
// 			(*psi).set(size,wf);
// 		}

// 	}


// }


// return the super-bosonic coherent state

MPS
super_bosonic_coherent_state( MPS *psi_coherent, const SiteSet sites_coherent, const complex<double> alpha, const double s, const double c, double dt, double T)
{
	int L = length(*psi_coherent);
	coherent_state_site_j( psi_coherent, sites_coherent, L , 1, alpha );
	adiabatic_transformation_linear_protocol( psi_coherent, sites_coherent, s, c, dt, T);

	return *psi_coherent;
}


// return the super-bosonic coherent state
/*
MPS
super_bosonic_coherent_state_from_ground_states( MPS *psi_coherent, const SiteSet sites, const complex<double> alpha, const double s, const double c, double dt, double T, Args const& physical_args, Args const& numerical_args)
{
	int L = length(*psi_coherent);
	int cut_off_fock_space = dim(sites(1))-1;
	double simmetry_sector;

	// GETTING SYMMETRY SECTOR FROM THE FILE

	ifstream symmetry_sector_file;
	string directory_symmetry = "/home/ricval/Documenti/Bosonic/Bosonic_Quantum_East_Model_Cpp/mmGcbQEM/symmetry_sector_minus";
	string name_symmetry_sector_file = format("%s/symmetry_sector%d_maxcutoff30_s%.2f_c%.2f.dat",directory_symmetry,0,s,c);
	cerr << name_symmetry_sector_file << endl;
	symmetry_sector_file.open(name_symmetry_sector_file);

	for (int i = 1; i < cut_off_fock_space; i++)
	{
		float tmp;
		symmetry_sector_file >> tmp;
	}

	symmetry_sector_file >> simmetry_sector;

	// BUILDING THE SUPER-COHERENT STATE VIA LINEAR COMBINATIONS OF GROUND STATES


	MPS super_coherent_state = randomMPS(sites);
	initial_state_vacuum_state_correct_link( &super_coherent_state, sites, size );
	super_coherent_state *= weigth_coherent_state( alpha, 0);
	for(int n0 = 1; n0 <= cut_off_fock_space ; n0++)
	{
		MPO H = H_mmGcbQEM_n0_untouched( sites, size , n0, simmetry_sector , s,  c);
		MPS ground_state = randomMPS(sites);
		initial_state_n0_excitation( &ground_state, sites, size , n0, 1);
		int bond_dimension = perform_DMRG( &ground_state ,  H , sites , 15 , physical_args, numerical_args);
		ground_state /= norm(ground_state);
		ground_state *= weigth_coherent_state( alpha, n0);
		super_coherent_state = sum(super_coherent_state, ground_state,{"MaxDim",500,"Cutoff",1E-8});
	}

	*psi_coherent = super_coherent_state;
	return *psi_coherent;
}

*/



// Return a super-bosonic squeezed state
MPS
super_bosonic_squeezed_state( MPS *psi_squeezed, const SiteSet sites_squeezed, const double alpha, const double s, const double c, double dt, double T)
{
	int L = length(*psi_squeezed);
	squeezed_state_site_j( psi_squeezed, sites_squeezed, L , 1, alpha );
	adiabatic_transformation_linear_protocol( psi_squeezed, sites_squeezed, s, c, dt, T);

	return *psi_squeezed;

}

void
kink_state(MPS *psi, const SiteSet sites, const int number_ones)
{
	int L =  length(*psi);
	
	for(int j=1 ; j<= number_ones; j++)     put_occupation(psi,sites,j,1);
	for(int j=number_ones + 1 ; j<= L; j++) put_occupation(psi,sites,j,0);

}

void
put_occupation(MPS *psi, const SiteSet sites, const int position, const int n)
{
	Index sj;
	Index li;
	Index ri;
	ITensor wf;
	int L = length(*psi);
	if( position == 1)
	{
		sj = sites(1);
		li = commonIndex((*psi)(1),(*psi)(2));
		wf = ITensor(sj,li);
		for( int d=1; d <= n; d++) wf.set(sj(d),li(1), 0);
		wf.set(sj(n+1),li(1),1);
		for( int d=n+2; d <= dim(sj); d++) wf.set(sj(d),li(1), 0);
	}

	else if( position == L)
	{
		sj  = sites(L);
		li = commonIndex((*psi)(L-1),(*psi)(L));
		wf = ITensor(sj,li);

		for( int d=1; d <= n; d++) wf.set(sj(d),li(1), 0);
		wf.set(sj(n+1),li(1),1);
		for( int d=n+2; d <= dim(sj); d++) wf.set(sj(d),li(1), 0);
	}

	else
	{
		sj = sites(position);
		li = commonIndex((*psi)(position-1),(*psi)(position));
		ri = commonIndex((*psi)(position),(*psi)(position+1));
		wf = ITensor(sj,li,ri);

		for( int d=1; d <= n; d++) wf.set(sj(d),li(1),ri(1), 0);
		wf.set(sj(n+1),li(1),ri(1),1);
		for( int d=n+2; d <= dim(sj); d++) wf.set(sj(d),li(1),ri(1), 0);
		
	}

	(*psi).set(position,wf);

}


// factorial n!

double
factorial(int n)
{
	double factorial = 1;
	for(int j=1; j<=n ; j++) factorial *= j;
	return factorial;	
}
	
// weight coherent state

complex<double> 
weigth_coherent_state( const complex<double> alpha, const int k)
{
	complex<double> alpha_power_k = pow(alpha,k);
	return exp(-abs(alpha)*abs(alpha) / 2 ) * (alpha_power_k.real() + 1i*alpha_power_k.imag()) / sqrt(factorial(k));
}

// weight squeezed state

double 
weigth_squeezed_state( const double r, const int k)
{
	double result;
	result  = pow((-1 * tanh(r)),int(k/2));
	result *= sqrt(factorial(k));
	result /= (pow(2,k/2) * factorial(int(k/2)));
	result /= sqrt(cosh(r));
	return result;
}


MPO
bQEM_dressed_operator(MPO O,  const SiteSet sites, double s, double c, double dt, double T)
{

	int size = length(sites);
	double J_target = exp(-s);
	int number_step = int(T/dt);

	for(int step=1; step<=int(number_step); step++)
	{
		double J = J_target * step * dt / T;
		cerr << J << endl;
		// exp_H_mmGcbQEM_n0_notfixed( const SiteSet sites, int size , int n0, double symmetry , double J, double c, double dt)

		MPO expH  = exp_H_mmGcbQEM_n0_notfixed( sites, size , 1, 0.0 , J,  c, dt);
		MPO expHd = exp_H_mmGcbQEM_n0_notfixed( sites, size , 1, 0.0 , J,  c, -1*dt);
		
		O = nmultMPO( O , prime(expH) ,{"MaxDim",1000,"Cutoff",1E-14}); 
		O.mapPrime(2,1);
		O = nmultMPO( expHd , prime(O) ,{"MaxDim",1000,"Cutoff",1E-14}); 
		O.mapPrime(2,1);

	}	
	return O;
}

