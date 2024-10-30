#include "observables.h"
#include "build_hamiltonian.h"
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
entanglement_entropy( MPS* psi , int site)
	{
	(*psi).position(site); 
	ITensor wf = (*psi)(site) * (*psi)(site+1);
	ITensor U  = (*psi)(site);
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

// void 
// measure_mx_mz( const SpinHalf sites , MPS psi , const int N)
// 	{
	
// 	for( int j = 1 ; j <= N ; j++ )
// 		{
// 		psi.position(j);
// 		Real Mx1 = 2 * (dag(prime(psi.A(j), "Site" )) * sites.op("Sx",j) * psi.A(j)).real();
// 		Real Mz1 = 2 * (dag(prime(psi.A(j), "Site" )) * sites.op("Sz",j) * psi.A(j)).real();
// 		cout << "Sx_" << j << " = " << Mx1 << "\n"
// 			 << "Sz_" << j << " = " << Mz1 << endl;
// 		}
// 	}

//----------------------------------------------------------------------
// expectation value: <\sigma_j^x^2>
double 
expectation_value_sigma_x_square( MPS *state , const SiteSet sites , const int j )
{
	ITensor observable = op(sites, "A", j);
	observable += op(sites,"Adag",j);
	ITensor observable2 = prime(observable);

	(*state).position(j);
	
	auto ket = (*state)(j);
	auto bra = dag(prime(prime(ket,"Site"),"Site"));
	double expectation_value = elt(bra * observable * observable2 * ket);

	return expectation_value;

}

//----------------------------------------------------------------------
// expectation value: <\sigma_j^x>
double 
expectation_value_sigma_x( MPS *state , const SiteSet sites , const int j )
{
	ITensor observable = op(sites, "A", j);
	observable += op(sites,"Adag",j);

	(*state).position(j);
	
	auto ket = (*state)(j);
	auto bra = dag(prime(ket,"Site"));
	double expectation_value = elt(bra * observable * ket);

	return expectation_value;

}


//----------------------------------------------------------------------
// expectation value: <\sigma_j^x n_j>
double 
expectation_value_sigma_x_n( MPS *state , const SiteSet sites , const int j )
{
	ITensor observable = op(sites, "A", j);
	observable += op(sites,"Adag",j);
	ITensor observable2 = prime(op(sites,"N",j));

	(*state).position(j);
	
	auto ket = (*state)(j);
	auto bra = dag(prime(prime(ket,"Site"),"Site"));
	double expectation_value = elt(bra * observable * observable2 * ket);

	return expectation_value;

}


//----------------------------------------------------------------------
// expectation value: <n_j \sigma_j^x>
double 
expectation_value_n_sigma_x( MPS *state , const SiteSet sites , const int j )
{
	ITensor observable = op(sites,"N",j);
	ITensor observable2 = op(sites, "A", j);
	observable2 += op(sites,"Adag",j);
	observable2 = prime(observable2);

	(*state).position(j);
	
	auto ket = (*state)(j);
	auto bra = dag(prime(prime(ket,"Site"),"Site"));
	double expectation_value = elt(bra * observable * observable2 * ket);

	return expectation_value;

}


// measure occupation number along the 1D chain

void 
measure_occupation_number( MPS *ground_state , const SiteSet sites , const int size ,  vector<double> &occupation_number )
{
for(int j = 1 ; j <= size ; j++)
	{
	ITensor observable = op(sites, "N" , j);
	(*ground_state).position(j);
	auto ket = (*ground_state)(j);
	auto bra = dag(prime(ket,"Site"));
	complex<double> expectation_value = eltC(bra * observable * ket);
	occupation_number.push_back( expectation_value.real() );
	}
}

double
measure_imbalance( vector<double> &occupation_number, int k)
{

	double nk = occupation_number[k];
	occupation_number.erase(occupation_number.begin()+k);
	double nmax = *max_element(occupation_number.begin(), occupation_number.end());
	
	return (nk-nmax)/(nk+nmax);

}


double max_projector_at_cutoff(vector<vector<double> > &projector_all_sites,const int size,const int cut_off)
{
	double maximum = -1;
	for(int j = 1 ; j <= size ; j++) maximum = std::max(projector_all_sites[j-1][cut_off-1], maximum);

	return maximum;

}

// measure squareoccupation number along the 1D chain

void 
measure_square_occupation_number( MPS *ground_state , const SiteSet sites , const int size ,  vector<double> &square_occupation_number )
{
for(int j = 1 ; j <= size ; j++)
	{
	ITensor observable = op(sites, "N" , j);
	ITensor observable2 = prime(op(sites, "N" , j));

	(*ground_state).position(j);
	auto ket = (*ground_state)(j);
	auto bra = dag(prime(prime(ket,"Site"),"Site"));
	double expectation_value = elt(bra * observable * observable2 * ket);
	square_occupation_number.push_back( expectation_value );
	}
	
}



//---------------------------------------------------------------------

// measure of the projector along all the sites and all the Fock space

void measure_projector_all_sites( MPS *ground_state , const SiteSet sites , const int size , const int cut_off_fock_space ,  vector<vector<double> > &projector_all_sites ,  vector<double> &occupation_number)
{
	for(int j = 1 ; j <= size ; j++)
	{
	vector<double> projector_single_site;
	(*ground_state).position(j);

	// cerr << "Norm ground state : " << norm(ground_state) << endl;

	double check_number_particles = 0.;
	double check_normalization_projector = 0. ;

	for( int n = 0 ; n <= cut_off_fock_space ; n++ )
		{
		Index physical_index = sites(j);
		Index physical_index_prime = prime(sites(j));
		ITensor projector = ITensor(physical_index, physical_index_prime); //initialized with all the elements equal to zero.
	
		projector.set(physical_index(n+1),physical_index_prime(n+1),1.);

		ITensor ket = (*ground_state)(j);
		ITensor bra = dag(prime((*ground_state)(j),"Site"));
		
		complex<double> expectation_value_c = eltC(bra * projector * ket);
		double expectation_value = expectation_value_c.real();
		projector_single_site.push_back(expectation_value);

		check_number_particles += n * expectation_value ;
		check_normalization_projector += expectation_value;
		//verifing projector
		// double occupation = elt(dag(prime(projector,"Site")) * op(sites, "N" , j) * projector);
		// cerr << "<N> (expected-computer): " << ( n - occupation ) << endl;
	
		}

	cerr << "Site : " << j << " Difference (occupation_number - projector) " << (occupation_number[j-1]-check_number_particles) << endl;	
	cerr << "Site : " << j << " Normalization projector (sum of projectors should be 1) : " << check_normalization_projector << endl;	
	projector_all_sites.push_back(projector_single_site);
	}
}

// Measure the covariance matrix size X size. Each element is <N_i N_j>_c - 
/*
void
measure_covariance_matrix_number_operator( MPS *psi , const SiteSet sites , vector<vector<double> > &covariance_matrix_NN)
{
	int size = length(*psi);

	(*psi).position(1);
	(*psi).normalize();
	
	for(int i = 1; i < 2; i++)
	{
		for(int j = i; j < i+1 ; j++)
		{

			ITensor Ni_op  = op(sites, "N" , i);
			ITensor Ai_op  = op(sites, "A" , i);
			ITensor Adi_op = op(sites, "Adag" , i);
			
			ITensor Nj_op  = op(sites, "N" , j);
			ITensor Aj_op  = op(sites, "A" , j);
			ITensor Adj_op = op(sites, "Adag" , j);


					ITensor adag2 = op(sites, "Adag", j) * prime(op(sites, "Adag",j));
		adag2.mapPrime(2,1);


			double NiNj;
			complex<double> Aid_Ajd ;
			complex<double> Aid_Aj  ;
			complex<double> Ai_Aj   ; 
			complex<double> Ai_Ajd  ;

			if( i != j)
			{
				NiNj    = abs(compute_two_point(psi,sites ,Ni_op,Nj_op , i , j));
				Aid_Ajd = compute_two_point(psi,sites,Adi_op,Adj_op , i , j);
				Ai_Aj   = compute_two_point(psi,sites,Ai_op,Aj_op , i , j);
				Aid_Aj  = compute_two_point(psi,sites,Adi_op,Aj_op  , i , j);
				Ai_Ajd  = compute_two_point(psi,sites,Ai_op,Adj_op  , i , j);	
			}


			else
			{
				cerr << "sono qui" << endl;	
				(*psi).position(i);
				auto ket = (*psi)(i);
				auto bra =  dag(prime(prime(ket,"Site"),"Site"));
				NiNj = elt(bra * Ni_op * prime(Ni_op) * ket);
				Aid_Ajd = eltC(bra * Adi_op * prime(Adi_op) * ket);
				Ai_Aj   = eltC(bra * Ai_op * prime(Ai_op) * ket);
				Aid_Aj = eltC(bra * Adi_op * prime(Ai_op) * ket);
				Ai_Ajd = eltC(bra * Ai_op * prime(Adi_op) * ket);

			}
		
			(*psi).position(i);
			auto ket = (*psi)(i);
			auto bra = dag(prime(ket,"Site"));
			double Ni = elt( bra * Ni_op * ket);
			complex<double> Ai = eltC( bra * Ai_op * ket);

			(*psi).position(j);
			ket = (*psi)(j);
			bra = dag(prime(ket,"Site"));
			double Nj = elt( bra * Nj_op * ket);
			complex<double> Aj = eltC( bra * Aj_op * ket);
		
			cerr << " Observables " << endl;
			cerr << "NN : " << NiNj << endl;
			cerr << "N : " << Ni << endl;
			cerr << "A : " << Ai << endl; 
			cerr << "AdA : " << Aid_Aj << endl;
			cerr << "AdAd : " << Aid_Ajd << endl;

//  CI STA UN PROBLEMA CON LA PRESENZA DI A_I AD_J . FORSE L'ORDINE E' SBAGLIATO. SULLO STATO DI VUOTO DEVO AVERE 0. DATO CHE E' AUTOSTATO DI N
	// CON STATI SQUEEZED VIENE VARIANZA NEGATIVA. IO DIREI DI TAGLIARE LA TESTA LA TORO 
	// SHIFTANDO GLI OPERATORI DI ANNICHILAZIONE E DISTRUZIONE DEL LORO VALORE ATTESO. COSI' WICK E' BANALE 
			complex<double> NiNj_gaussian_approx_connected = Aid_Ajd * Ai_Aj + Aid_Aj * Ai_Aj - 2 * abs(Ai) * abs(Ai) * abs(Aj) * abs(Aj);
			double NiNj_connected = NiNj - Ni*Nj;
			double coherence = abs(Ai)*abs(Ai) - Ni;
			// cerr << "coherence " << i << " " << coherence << endl; 
			cerr << "N_" << i << "N_" << j << " : " << NiNj_gaussian_approx_connected << " " << NiNj_connected << endl;
		}
	}
}


*/
void
measure_covariance_matrix_number_operator( MPS *psi , const SiteSet sites , vector<vector<double> > &covariance_matrix_NN_system, vector<vector<double> > &covariance_matrix_NN_gaussian, vector<vector<double> > &relative_error)
{
	int size = length(*psi);

	(*psi).position(1);
	(*psi).normalize();
	
	for(int i = 1; i <= size; i++)
	{
		vector<double> covariance_ij;
		vector<double> covariance_ij_gaussian;
		vector<double> relative_error_ij;

		for(int j = 1; j <= size ; j++)
		{
			
			double NiNj;
			complex<double> Aid_Ajd ;
			complex<double> Aid_Aj  ;
			complex<double> Ai_Aj   ; 
			complex<double> Ai_Ajd  ;

			// define A 
			ITensor Ai_op  = op(sites, "A" , i);
			ITensor Aj_op  = op(sites, "A" , j);
			ITensor Adi_op  = op(sites, "Adag" , i);
			ITensor Adj_op  = op(sites, "Adag" , j);

			(*psi).position(i);
			auto ket = (*psi)(i);
			auto bra = dag(prime(ket,"Site"));
			complex<double> Ai = eltC( bra * Ai_op * ket);

			(*psi).position(j);
			ket = (*psi)(j);
			bra = dag(prime(ket,"Site"));
			// bra = prime(ket,"Site");
			// bra = prime(bra,"Site");
			// bra = dag(bra);
			complex<double> Aj = eltC( bra * Aj_op * ket);

			// attempt - 14.12.21
			ITensor Ni_op  = op(sites, "N" , i);
			ITensor Nj_op  = op(sites, "N" , j);

			// end attempt - 14.12.21
			// Adi_op -= conj(Ai) * op(sites,"Id",i);
			// Ai_op  -= Ai * op(sites,"Id",i);

			// Adj_op -= conj(Aj) * op(sites,"Id",j);
			// Aj_op  -= Aj * op(sites,"Id",j);

			// commented from here - 14.12.21
			// ITensor Ni_op  = Adi_op * prime(Ai_op);
			// ITensor Nj_op  = Adj_op * prime(Aj_op);

			// Ni_op.mapPrime(2,1);
			// Nj_op.mapPrime(2,1);
			// to here - 14.12.21
			
			// PrintData(Ni_op);
			// exit(0);

			if( i != j)
			{
				NiNj    = abs(compute_two_point(psi,sites ,Ni_op,Nj_op , i , j));
				Aid_Ajd = compute_two_point(psi,sites,Adi_op,Adj_op , i , j);
				Ai_Aj   = compute_two_point(psi,sites,Ai_op,Aj_op , i , j);
				Aid_Aj  = compute_two_point(psi,sites,Adi_op,Aj_op  , i , j);
				Ai_Ajd  = compute_two_point(psi,sites,Ai_op,Adj_op  , i , j);	
			}


			else
			{
				(*psi).position(i);
				auto ket = (*psi)(i);
				auto bra =  dag(prime(prime(ket,"Site"),"Site"));
				NiNj = elt(bra * prime(Ni_op) * Ni_op * ket);
				Aid_Ajd = eltC(bra * prime(Adi_op) * Adi_op * ket);
				Ai_Aj   = eltC(bra * prime(Ai_op) * Ai_op * ket);
				Aid_Aj = eltC(bra * prime(Adi_op) * Ai_op * ket);
				Ai_Ajd = eltC(bra * prime(Ai_op) * Adi_op * ket);

			}
		
			(*psi).position(i);
			ket = (*psi)(i);
			bra = dag(prime(ket,"Site"));
			double Ni = elt( bra * Ni_op * ket);
			Ai = eltC( bra * Ai_op * ket);
 
			(*psi).position(j);
			ket = (*psi)(j);
			bra = dag(prime(ket,"Site"));
			double Nj = elt( bra * Nj_op * ket);
			Aj = eltC( bra * Aj_op * ket);

			if( i == j && i==1)
			{
			cerr << " Observables " << endl;
			cerr << "NN : " << NiNj << endl;
			cerr << "N : " << Ni << endl;
			cerr << "A : " << Ai << endl; 
			cerr << "AdA : " << Aid_Aj << endl;
			cerr << "AdAd : " << Aid_Ajd << endl;
			}
			// complex<double> NiNj_gaussian_approx_connected = Aid_Ajd * Ai_Aj + Aid_Aj * Ai_Ajd ;
			// double NiNj_connected = NiNj - Ni*Nj;
			// 14.12.21 post Girvin: divide the covariance by < n_i n_j>.
			// In this way, along the diagonal we have 1 if the state is coherent (the variance is N).


			// measure average occupation number
			ITensor observable = op(sites, "N" , i);
			(*psi).position(i);
			ket = (*psi)(i);
			bra = dag(prime(ket,"Site"));
			double ni = eltC(bra * observable * ket).real();

			observable = op(sites, "N" , j);
			(*psi).position(j);
			ket = (*psi)(j);
			bra = dag(prime(ket,"Site"));
			double nj = eltC(bra * observable * ket).real();

			if( i == j && i==1)
			{
			cerr << "N_" << i << " " << ni << " " << Ni << endl;

			
			}




			complex<double> NiNj_gaussian_approx_connected = (Aid_Ajd * Ai_Aj + Aid_Aj * Ai_Ajd) - 2 * abs(Ai)*abs(Ai) * abs(Aj)*abs(Aj);
			double NiNj_connected = (NiNj - Ni*Nj);
		
			// NiNj_gaussian_approx_connected /= Ni;
			// NiNj_connected /=  Ni;	
			// if you want to have not-normalized things, simply comment the above two lines
			if( i == j && i==1)
			{
			cerr << "N_" << i << "N_" << j << " : " << NiNj_gaussian_approx_connected << " " << NiNj_connected << " " << (abs(NiNj_gaussian_approx_connected)-NiNj_connected)/abs(NiNj_gaussian_approx_connected) << endl;
			}
			covariance_ij.push_back(NiNj_connected);
			covariance_ij_gaussian.push_back(NiNj_gaussian_approx_connected.real());
			relative_error_ij.push_back((NiNj_gaussian_approx_connected.real()-NiNj_connected)/NiNj_connected);

		}
		covariance_matrix_NN_system.push_back(covariance_ij);
		covariance_matrix_NN_gaussian.push_back(covariance_ij_gaussian);
		relative_error.push_back(relative_error_ij);
		
	}
}

complex<double>
compute_two_point( MPS *psi, const SiteSet sites, ITensor op_i, ITensor op_j, int i, int j)
{
	if(j<i)
	{
		int k = i;
		ITensor op_k = op_i;
		i = j;
		j = k;
		op_i = op_j;
		op_j = op_k;
	}

	//'gauge' the MPS to site i
	//any 'position' between i and j, inclusive, would work here
	(*psi).position(i); 

	//Create the bra/dual version of the MPS psi
	auto psidag = dag(*psi);

	//Prime the link indices to make them distinct from
	//the original ket links
	psidag.prime("Link");

	//index linking i-1 to i:
	auto li_1 = leftLinkIndex(*psi,i);

	auto C = prime((*psi)(i),li_1)*op_i;
	C *= prime(psidag(i),"Site");
	for(int k = i+1; k < j; ++k)
		{
		C *= (*psi)(k);
		C *= psidag(k);
		}
	//index linking j to j+1:
	auto lj = rightLinkIndex((*psi),j);

	C *= prime((*psi)(j),lj)*op_j;
	C *= prime(psidag(j),"Site");

	complex<double> result = eltC(C); //or eltC(C) if expecting complex	
	return result;
}


double 
compute_variance_H_mmGcbQEM(MPS *psi , const SiteSet sites, int size , int cut_off_fock_space, int n0, int symmetry_sector, double s, double c)
{
	double symmetry;
	ifstream symmetry_sector_file;
	string directory_symmetry = "/home/ricval/Documenti/Bosonic/Bosonic_Quantum_East_Model_Cpp/mmGcbQEM/symmetry_sector_minus";
	string name_symmetry_sector_file = format("%s/symmetry_sector%d_maxcutoff30_s%.2f_c%.2f.dat",directory_symmetry,int(symmetry_sector),s,c);
	cerr << name_symmetry_sector_file << endl;
	symmetry_sector_file.open(name_symmetry_sector_file);

	for (int i = 1; i < cut_off_fock_space; i++)
	{
		float tmp;
		symmetry_sector_file >> tmp;
	}

	symmetry_sector_file >> symmetry;

	MPO H = H_mmGcbQEM( sites, size , n0, symmetry , s, c) ; 


	double variance = inner((*psi),H,H,(*psi))  -   inner((*psi),H,(*psi)) *  inner((*psi),H,(*psi));
	cerr << "I have computed variance " << endl;
	return variance;
}




double
measure_delta_x(MPS *psi, MPO *A, MPO *Adag )
{

	complex<double> a_expectation_value     = innerC(*psi, *A   , *psi);
	complex<double> aa_expectation_value    = innerC(*psi, *A   , *A, *psi);
	complex<double> adaga_expectation_value = innerC(*psi, *Adag, *A, *psi);

	complex<double> delta_x = 1 + 2 * adaga_expectation_value + 2 * aa_expectation_value.real() - 4 * a_expectation_value.real() * a_expectation_value.real();

	return delta_x.real();
}

double
measure_delta_p(MPS *psi, MPO *A, MPO *Adag )
{

	complex<double> a_expectation_value     = innerC(*psi, *A   , *psi);
	complex<double> aa_expectation_value    = innerC(*psi, *A   , *A, *psi);
	complex<double> adaga_expectation_value = innerC(*psi, *Adag, *A, *psi);

	complex<double> delta_p = 1 + 2 * adaga_expectation_value - 2 * aa_expectation_value.real() - 4 * a_expectation_value.imag() * a_expectation_value.imag();

	return delta_p.real();

}


double
measure_squeezing(MPS *psi, const SiteSet sites, const int j)
{
	if( j <= length(*psi) )
	{
		ITensor Nj = op(sites, "N" , j);
		ITensor adag2 = op(sites, "Adag", j) * prime(op(sites, "Adag",j));
		adag2.mapPrime(2,1);
		

		(*psi).position(j);
		auto ket = (*psi)(j);
		auto bra = dag(prime(ket,"Site"));
		complex<double> Nj_exp = eltC(bra * Nj * ket);
		complex<double> adag2_exp = eltC(bra * adag2 * ket);
		double squeezing = 	1 + 2 * Nj_exp.real() - 2 * abs(adag2_exp);
		return squeezing;
	}
	else return -1;
}

double
measure_dressed_squeezing(MPS *psi, const SiteSet sites, MPO A, MPO N)
{
	complex<double> N_measured = innerC((*psi),N,(*psi));
	complex<double> A2_measured = innerC((*psi),A,A,(*psi));

	double squeezing_dressed = 1 + 2 * N_measured.real() - 2 * abs(A2_measured);
	return squeezing_dressed;
}
