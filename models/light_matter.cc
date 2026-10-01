/**
 * @file light_matter.cc
 * @brief Implementation of light_matter.h (the functions are documented in the header).
 */
#include "light_matter.h"
#include <itensor/all.h>
#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

using namespace std;
using namespace itensor;


// -----------------------------------------------------------------
// Gates of the Tavis-Cummings model 
// H = omegao n_a + h S^z + g (S^+ a + S^- a^\dag)
// where we have collective operators S^alpha = \sum_{j=1}^{N-1} s_j^alpha 
// We consider (N-1) spin-1/2 explicitly (we do not use collective operators).

// It returns either short-range gates or long-range one, where
// H_short = omegao n_a + h S^z 
// H_long  = g (S^+ a + S^- a^\dag)
// We DO NOT implement the swap gates as gates, but directly in the TEBD algorithm.

vector<BondGate>
gates_tavis_cummings(const SiteSet sites , const double omega0 , const double h , const double g, const double dt, string matter_or_photon)
{
    vector<BondGate> gates;

	int N = length(sites);
	
	// -----------------------------------------------------------------------------------
	// gates for H = \omega_0 n_a + \sum_j h \sigma_j^z
	// single-site gate (it implies no need of SVD when applied)
	// in order to use the BondGate class I write as two-gate Hamiltonian
	// When applied, I can exploit the fact there is no need of SVDs

	if(matter_or_photon == "short-range")
	{
	cerr << "Build short-range interactions" << endl;

	for(int j=1 ; j<N;j++)
    {
        
        Index sj = sites(j);
        Index sjp = prime(sites(j));

		Index sj1 = sites(j+1);
        Index sj1p = prime(sites(j+1));
        
        ITensor S_j  = ITensor(sj ,sjp );
		ITensor I_j  = ITensor(sj ,sjp );
		
        ITensor S_j1 = ITensor(sj1,sj1p);
        ITensor I_j1 = ITensor(sj1,sj1p);

		vector<double> local_fields;
	
		if(hasTags(sj,"Site,Boson"))
		{
			for(int d=1; d <= dim(sj) ; d++) S_j.set(sj(d),sjp(d),d-1.);
			local_fields.push_back(omega0);
        }
				
		if(hasTags(sj,"Site,S=1/2"))
		{
			S_j.set(sj(1),sjp(1),1.);
			S_j.set(sj(2),sjp(2),-1.);	
			local_fields.push_back(h);
		}

		if(hasTags(sj1,"Site,Boson"))
		{
			for(int d=1; d <= dim(sj1) ; d++) S_j1.set(sj1(d),sj1p(d),d-1.);
			local_fields.push_back(omega0);

        }
				
		if(hasTags(sj1,"Site,S=1/2"))
		{
			S_j1.set(sj1(1),sj1p(1),1.);
			S_j1.set(sj1(2),sj1p(2),-1.);	
			local_fields.push_back(h);

		}

		for(int d=1; d <= dim(sj)  ; d++) I_j.set(sj(d),sjp(d),1.);
		for(int d=1; d <= dim(sj1) ; d++) I_j1.set(sj1(d),sj1p(d),1.);

		ITensor hj;
		if(j==1) 		hj = local_fields[0]    * S_j * I_j1 + local_fields[1]/2. * I_j * S_j1;
		else if(j==N-1) hj = local_fields[0]/2. * S_j * I_j1 + local_fields[1]    * I_j * S_j1;
		else 			hj = local_fields[0]/2. * S_j * I_j1 + local_fields[1]/2. * I_j * S_j1;
		BondGate g = BondGate(sites,j,j+1,BondGate::tReal,dt/2.,hj); 
		gates.push_back(g);

    }

	for(int j=N-1 ; j>=1;j--)
    {
        
        Index sj = sites(j);
        Index sjp = prime(sites(j));

		Index sj1 = sites(j+1);
        Index sj1p = prime(sites(j+1));
        
        ITensor S_j  = ITensor(sj ,sjp );
		ITensor I_j  = ITensor(sj ,sjp );
		
        ITensor S_j1 = ITensor(sj1,sj1p);
        ITensor I_j1 = ITensor(sj1,sj1p);

		vector<double> local_fields;
	
		if(hasTags(sj,"Site,Boson"))
		{
			for(int d=1; d <= dim(sj) ; d++) S_j.set(sj(d),sjp(d),d-1.);
			local_fields.push_back(omega0);
        }
				
		if(hasTags(sj,"Site,S=1/2"))
		{
			S_j.set(sj(1),sjp(1),1.);
			S_j.set(sj(2),sjp(2),-1.);	
			local_fields.push_back(h);
		}

		if(hasTags(sj1,"Site,Boson"))
		{
			for(int d=1; d <= dim(sj1) ; d++) S_j1.set(sj1(d),sj1p(d),d-1.);
			local_fields.push_back(omega0);

        }
				
		if(hasTags(sj1,"Site,S=1/2"))
		{
			S_j1.set(sj1(1),sj1p(1),1.);
			S_j1.set(sj1(2),sj1p(2),-1.);	
			local_fields.push_back(h);

		}

		for(int d=1; d <= dim(sj)  ; d++) I_j.set(sj(d),sjp(d),1.);
		for(int d=1; d <= dim(sj1) ; d++) I_j1.set(sj1(d),sj1p(d),1.);
		
		ITensor hj;
		if(j==1) hj = local_fields[0] * S_j * I_j1 + local_fields[1]/2. * I_j * S_j1;
		else if(j==N-1) hj = local_fields[0]/2. * S_j * I_j1 + local_fields[1] * I_j * S_j1;
		else hj = local_fields[0]/2. * S_j * I_j1 + local_fields[1]/2. * I_j * S_j1;
		BondGate g = BondGate(sites,j,j+1,BondGate::tReal,dt/2.,hj); 
		gates.push_back(g);

    }

	return gates;
	}
	// -----------------------------------------------------------------------------------

	// gate H = g \sum_j s_j^+ a + h.c.
	// I implement swap gate as gates using the built-in gate  BondGate(sites,7,8); (see https://itensor.org/docs.cgi?vers=cppv3&page=classes/bondgate)
	// NO: IT DOES NOT WORK WITH DIFFERENT PHYSICAL DIMENSIONS AS IT IS THE CASE OF SPIN-BOSON MODEL
	// Alternative: I implement the swap gate using SVDs -> drawback: I have to explicitly write it in the main

	else if(matter_or_photon == "long-range")
	{

	Index sph   = sites(1);
	Index sph_p = prime(sites(1));

	ITensor A  = ITensor(sph,sph_p);
	ITensor Ad = ITensor(sph,sph_p);
	
	for(int d=1; d < dim(sph) ; d++)
	{
		A.set( sph(d)  ,sph_p(d+1),sqrt(d));
		Ad.set(sph(d+1),sph_p(d)  ,sqrt(d));
	}
	cerr << "Build photon-matter terms" << endl;

	for(int j=1 ; j<N;j++)
    {
		Index sj1   = sites(j+1);
        Index sj1p  = prime(sites(j+1));
        ITensor Sp = ITensor(sj1,sj1p);
        ITensor Sm = ITensor(sj1,sj1p);


		if(hasTags(sj1,"Site,S=1/2"))
		{
			Sm.set(sj1(1),sj1p(2),1.);
			Sp.set(sj1(2),sj1p(1),1.);	
		}


		ITensor hj = g * (A * Sp + Ad * Sm);
		BondGate g = BondGate(sites,1,j+1,BondGate::tReal,dt/2.,hj); 
		gates.push_back(g);

    }

	for(int j=N-1 ; j>=1;j--)
    {
		Index sj1   = sites(j+1);
        Index sj1p  = prime(sites(j+1));
        ITensor Sp = ITensor(sj1,sj1p);
        ITensor Sm = ITensor(sj1,sj1p);


		if(hasTags(sj1,"Site,S=1/2"))
		{
			Sm.set(sj1(1),sj1p(2),1.);
			Sp.set(sj1(2),sj1p(1),1.);	
		}


		ITensor hj = g * (A * Sp + Ad * Sm);
		BondGate g = BondGate(sites,1,j+1,BondGate::tReal,dt/2.,hj); 
		gates.push_back(g);
    }
	return gates;
	}

	else{
		cerr << "string argument inserted is neither 'short-range' nor 'long-range" << endl;
		cerr << "Return empy gates" << endl;
		return gates;
	}
}


// -----------------------------------------------------------------
// Gates of the general light-matter interaction. It is possible to select:
// type_of_coupling = "dicke" : dicke type interactions           -> g S^x (a+a^\dag)
// type_of_coupling = "tavis" : jaynes cummings type interactions -> g (S^+ a + S^- a^\dag)

// It returns either short-range gates or long-range one, where
// H_short = omegao n_a + h S^z 
// H_long  = dicke or tavis
// We DO NOT implement the swap gates as gates, but directly in the TEBD algorithm.

vector<BondGate>
gates_photon_matter(const SiteSet sites , const double omega0 , const double h , const double g, const double dt, string matter_or_photon, string type_of_coupling, const double V, string interaction_axis)
{
    vector<BondGate> gates;

	int N = length(sites);
	
	// -----------------------------------------------------------------------------------
	// gates for H = \omega_0 n_a + \sum_j h \sigma_j^z
	// single-site gate (it implies no need of SVD when applied)
	// in order to use the BondGate class I write as two-gate Hamiltonian
	// When applied, I can exploit the fact there is no need of SVDs

	if(matter_or_photon == "short-range")
	{
		cerr << "Build short-range interactions" << endl;

		for(int j=1 ; j<N;j++)
		{
			
			Index sj = sites(j);
			Index sjp = prime(sites(j));

			Index sj1 = sites(j+1);
			Index sj1p = prime(sites(j+1));
			
			ITensor Sx_j = ITensor(sj ,sjp );
			ITensor Sz_j = ITensor(sj ,sjp );
			ITensor I_j  = ITensor(sj ,sjp );
			
			ITensor Sx_j1 = ITensor(sj1,sj1p);
			ITensor Sz_j1 = ITensor(sj1,sj1p);
			ITensor I_j1  = ITensor(sj1,sj1p);

			vector<double> local_fields;
		
			if(hasTags(sj,"Site,Boson"))
			{
				for(int d=1; d <= dim(sj) ; d++) Sz_j.set(sj(d),sjp(d),d-1.);
				local_fields.push_back(omega0);
			}
					
			if(hasTags(sj,"Site,S=1/2"))
			{
				Sz_j.set(sj(1),sjp(1),1.);
				Sz_j.set(sj(2),sjp(2),-1.);	
				Sx_j.set(sj(1),sjp(2),1.);
				Sx_j.set(sj(2),sjp(1),1.);	
				local_fields.push_back(h);
			}

			if(hasTags(sj1,"Site,Boson"))
			{
				for(int d=1; d <= dim(sj1) ; d++) Sz_j1.set(sj1(d),sj1p(d),d-1.);
				local_fields.push_back(omega0);

			}
					
			if(hasTags(sj1,"Site,S=1/2"))
			{
				Sz_j1.set(sj1(1),sj1p(1),1.);
				Sz_j1.set(sj1(2),sj1p(2),-1.);	
				Sx_j1.set(sj1(1),sj1p(2),1.);
				Sx_j1.set(sj1(2),sj1p(1),1.);	
				local_fields.push_back(h);

			}

			for(int d=1; d <= dim(sj)  ; d++) I_j.set(sj(d),sjp(d),1.);
			for(int d=1; d <= dim(sj1) ; d++) I_j1.set(sj1(d),sj1p(d),1.);

			// local fields diagonal in Sz and Nphoton
			ITensor hj;
			if(j==1) 		hj = local_fields[0]    * Sz_j * I_j1 + local_fields[1]/2. * I_j * Sz_j1;
			else if(j==N-1) hj = local_fields[0]/2. * Sz_j * I_j1 + local_fields[1]    * I_j * Sz_j1;
			else 			hj = local_fields[0]/2. * Sz_j * I_j1 + local_fields[1]/2. * I_j * Sz_j1;

			// if the sites considered both contain spin degrees of freedom
			if( hasTags(sj,"Site,S=1/2") &&  hasTags(sj1,"Site,S=1/2"))
			{
				ITensor h_nn ; 
				if(interaction_axis == "z")
				{
					cerr << "Building NN interactions.\n";
					ITensor Nj  = (I_j  - Sz_j )/2.;
					ITensor Nj1 = (I_j1 - Sz_j1)/2.;
					h_nn = V * Nj * Nj1;
				}
				else if(interaction_axis == "x")
				{
					cerr << "Building XX interactions.\n";
					h_nn = V * Sx_j * Sx_j1;
				}
				else
				{
					cerr << "Short-range interactions along " << interaction_axis << " not yet implemented.\n";
					exit(0);
				}
				hj += h_nn;
			}

			BondGate g = BondGate(sites,j,j+1,BondGate::tReal,dt/2.,hj); 
			gates.push_back(g);

		}

		for(int j=N-1 ; j>=1;j--)
		{
			
			Index sj = sites(j);
			Index sjp = prime(sites(j));

			Index sj1 = sites(j+1);
			Index sj1p = prime(sites(j+1));
			
			ITensor Sx_j = ITensor(sj ,sjp );
			ITensor Sz_j = ITensor(sj ,sjp );
			ITensor I_j  = ITensor(sj ,sjp );
			
			ITensor Sx_j1 = ITensor(sj1,sj1p);
			ITensor Sz_j1 = ITensor(sj1,sj1p);
			ITensor I_j1  = ITensor(sj1,sj1p);

			vector<double> local_fields;
		
			if(hasTags(sj,"Site,Boson"))
			{
				for(int d=1; d <= dim(sj) ; d++) Sz_j.set(sj(d),sjp(d),d-1.);
				local_fields.push_back(omega0);
			}
					
			if(hasTags(sj,"Site,S=1/2"))
			{
				Sz_j.set(sj(1),sjp(1),1.);
				Sz_j.set(sj(2),sjp(2),-1.);	
				Sx_j.set(sj(1),sjp(2),1.);
				Sx_j.set(sj(2),sjp(1),1.);	
				local_fields.push_back(h);
			}

			if(hasTags(sj1,"Site,Boson"))
			{
				for(int d=1; d <= dim(sj1) ; d++) Sz_j1.set(sj1(d),sj1p(d),d-1.);
				local_fields.push_back(omega0);

			}
					
			if(hasTags(sj1,"Site,S=1/2"))
			{
				Sz_j1.set(sj1(1),sj1p(1),1.);
				Sz_j1.set(sj1(2),sj1p(2),-1.);
				Sx_j1.set(sj1(1),sj1p(2),1.);
				Sx_j1.set(sj1(2),sj1p(1),1.);	
				local_fields.push_back(h);

			}

			for(int d=1; d <= dim(sj)  ; d++) I_j.set(sj(d),sjp(d),1.);
			for(int d=1; d <= dim(sj1) ; d++) I_j1.set(sj1(d),sj1p(d),1.);
			
			ITensor hj;
			if(j==1) 		hj = local_fields[0]    * Sz_j * I_j1 + local_fields[1]/2. * I_j * Sz_j1;
			else if(j==N-1) hj = local_fields[0]/2. * Sz_j * I_j1 + local_fields[1]    * I_j * Sz_j1;
			else 			hj = local_fields[0]/2. * Sz_j * I_j1 + local_fields[1]/2. * I_j * Sz_j1;

			// if the sites considered both contain spin degrees of freedom
			if( hasTags(sj,"Site,S=1/2") &&  hasTags(sj1,"Site,S=1/2"))
			{
				ITensor h_nn ; 
				if(interaction_axis == "z")
				{
					cerr << "Building NN interactions.\n";
					ITensor Nj  = (I_j  - Sz_j )/2.;
					ITensor Nj1 = (I_j1 - Sz_j1)/2.;
					h_nn = V * Nj * Nj1;
				}
				else if(interaction_axis == "x")
				{
					cerr << "Building XX interactions.\n";
					h_nn = V * Sx_j * Sx_j1;
				}	
				else
				{
					cerr << "Short-range interactions along " << interaction_axis << " not yet implemented.\n";
					exit(0);
				}
				hj += h_nn;
			}

			BondGate g = BondGate(sites,j,j+1,BondGate::tReal,dt/2.,hj); 
			gates.push_back(g);

		}

		return gates;
	
	}
	// -----------------------------------------------------------------------------------

	// gate H = g \sum_j s_j^+ a + h.c.
	// I implement swap gate as gates using the built-in gate  BondGate(sites,7,8); (see https://itensor.org/docs.cgi?vers=cppv3&page=classes/bondgate)
	// NO: IT DOES NOT WORK WITH DIFFERENT PHYSICAL DIMENSIONS AS IT IS THE CASE OF SPIN-BOSON MODEL
	// Alternative: I implement the swap gate using SVDs -> drawback: I have to explicitly write it in the main
	// photon operators

	else if(matter_or_photon == "long-range")
	{

	Index sph   = sites(1);
	Index sph_p = prime(sites(1));

	ITensor A  = ITensor(sph,sph_p);
	ITensor Ad = ITensor(sph,sph_p);
	
	for(int d=1; d < dim(sph) ; d++)
	{
		Ad.set( sph(d)  ,sph_p(d+1),sqrt(d));
		A.set(sph(d+1),sph_p(d)  ,sqrt(d));
	}
	cerr << "Build photon-matter terms" << endl;

	for(int j=1 ; j<N;j++)
    {
		Index sj1   = sites(j+1);
        Index sj1p  = prime(sites(j+1));
        ITensor Sp = ITensor(sj1,sj1p);
        ITensor Sm = ITensor(sj1,sj1p);


		if(hasTags(sj1,"Site,S=1/2"))
		{
			Sm.set(sj1(1),sj1p(2),1.);
			Sp.set(sj1(2),sj1p(1),1.);	
		}

		ITensor hj;
		if(type_of_coupling == "dicke")       hj = g * (A + Ad) * (Sp + Sm);
		else if (type_of_coupling == "tavis") hj = g * (A * Sp + Ad * Sm);
		else
		{
			cerr << "Selected type_of_coupling : " << type_of_coupling << " not supported.\n";
			exit(-1); 
		}
		BondGate g = BondGate(sites,1,j+1,BondGate::tReal,dt/2.,hj); 
		gates.push_back(g);

    }

	for(int j=N-1 ; j>=1;j--)
    {
		Index sj1   = sites(j+1);
        Index sj1p  = prime(sites(j+1));
        ITensor Sp = ITensor(sj1,sj1p);
        ITensor Sm = ITensor(sj1,sj1p);


		if(hasTags(sj1,"Site,S=1/2"))
		{
			Sm.set(sj1(1),sj1p(2),1.);
			Sp.set(sj1(2),sj1p(1),1.);	
		}


		ITensor hj;
		if(type_of_coupling == "dicke")       hj = g * (A + Ad) * (Sp + Sm);
		else if (type_of_coupling == "tavis") hj = g * (A * Sp + Ad * Sm);
		else
		{
			cerr << "Selected type_of_coupling : " << type_of_coupling << " not supported.\n";
			exit(-1); 
		}
		BondGate g = BondGate(sites,1,j+1,BondGate::tReal,dt/2.,hj); 
		gates.push_back(g);
    }
	return gates;
	}

	else{
		cerr << "string argument inserted is neither 'short-range' nor 'long-range" << endl;
		cerr << "Return empy gates" << endl;
		return gates;
	}
}
