//Note: sigma_i = 2 S_i

#include <itensor/all.h>
#include <sys/stat.h>
#include <iostream>
#include <string>
#include <vector>
#include <ctime>
#include <fstream>	//output file
#include <sstream>	//for ostringstream
#include <iomanip>
#include "mytn.h"
#include <filesystem>

using namespace std;
using namespace itensor;
namespace fs = std::filesystem;

//------------------------------------------------------------------MAIN--------------------------------------------------------------------

int main(int argc , char* argv[]){
		
	//getting input
	if(argc != 11)
		{
		cout << "Input: " << argv[0] << "<INITIAL STATE, UP(0), DOWN(1), DOMAIN-WALL(2)> <N> <J> <hxChoice> <hzChoice> <ttotal> <tstep> <nmeas> <bonddim> <LOCALE (0) / CLUSTER (1)>" << endl;
		return -1;
		}
	int state , N, nmeas, bonddim , localvscluster;
	int hxChoice = atoi( argv[3] );
	int hzChoice = atoi( argv[4] );	
	double J , hx , hz , ttotal , tstep ;
	double wall_time = 12. * 60 * 60;  //hours * minutes * seconds
		
	// all files are read/written in the folder data/
	fs::create_directories("data");
	fs::current_path("data");

	get_data( argv , &state , &N , &J , &hx , &hz , &ttotal, &tstep , &nmeas , &bonddim, &localvscluster);
	print_input( N , J , hx , hz , ttotal , tstep , nmeas , bonddim );
	
	cout << "Initial state ( UP = 0 , DOWN = 1 , DOMAIN WALL = 2) : " << state << endl;
	//Defining Sites and random MPS
	SpinHalf sites = SpinHalf(N,{"ConserveQNs=",false});
		
	//creating initial state: product states polarized along x
	string config;
	if( state == 0 )      config = string(N,'0');                          // all |+x>
	else if( state == 1 ) config = string(N,'1');                          // all |-x>
	else if( state == 2 ) config = string(N/2,'0') + string(N - N/2,'1');  // domain wall
	else { cerr << "Unknown initial state " << state << endl; return -1; }
	MPS psi = initial_computational_state( sites , config , "x" );
	
	//testing some observables in order to check the initial state.	
	measure_mx_mz( sites , psi , N);
	
	//TEBD

	//parameters defining the truncation during time evolution
	auto args = Args("Cutoff=",1E-16,"Verbose=",false , "MaxDim=" , bonddim );		
	
	int numberStep = int( ttotal/tstep + (1E-9*(ttotal/tstep)) ); 
	
	//building quenched Hamiltonian
	
	//Define the type "Gate" as a shorthand for BondGate<ITensor>
	using Gate = BondGate;   
	
	//Create a std::vector to hold the Trotter gates 
	auto gates = vector<Gate>();
	for(int b = 1; b <= N-1; b++)
		{
		ITensor hterm;
		build_single_step( &hterm , sites , N , J , hx , hz , b );
		Gate g = Gate(sites,b,b+1,Gate::tReal,tstep/2.,hterm); 
		gates.push_back(g);
		}
	for(int b = N-1; b >= 1; b-=1)
		{
		ITensor hterm;
		build_single_step( &hterm , sites , N , J , hx , hz , b );
		Gate g = Gate(sites,b,b+1,Gate::tReal,tstep/2.,hterm); 
		gates.push_back(g);
		}
		
	//file where save everything
	stringstream psi_file , sites_file , out;
	build_file_TEBD( &sites_file , &psi_file , N );

	//saving step 0 and input.txt
	out << psi_file.str() << "0";
	writeToFile(sites_file.str() , sites);
	writeToFile(out.str() , psi);	
	
	//starting time evolution
	cout << "Starting Real Time-Evolution using Trotter decomposition." << endl;
	cout << "Wall time : " << wall_time << endl;
	
	time_t c_begin = time(NULL);
	time_t c_step_start;
	time_t c_step_end;
	
	string entropy_file = "entropy.dat";
	ofstream SaveEntropy( entropy_file.c_str() , ios::app );
	
	for( int n = 1 ; n <= numberStep / nmeas ; n++)
		{
		c_step_start = time(NULL);
		gateTEvol( gates , nmeas * tstep , tstep , psi , args); //I want to save the state every nmeas steps.
		
		out.str("");
		out << psi_file.str() << n * nmeas; 
		cout << "Saving state at time" << n * nmeas * tstep << endl;
		writeToFile(out.str(), psi);
		cout << "State saved!" << endl;
		
		c_step_end = time(NULL);
		time_t time_elapsed_step = c_step_end - c_step_start;
		time_t time_elapsed_total = c_step_end - c_begin;
		double entropy = print_info( time_elapsed_step , time_elapsed_total , &psi , N , nmeas , n , tstep );
		SaveEntropy << n * nmeas * tstep  << " " << entropy << "\n";
		}

	cout << "Completed Real Time-Evolution using Trotter decomposition.\n"
		 << "You have simulated the following system.\n"
		 << "N: " << N << "\n"
		 << "J: " << J << "\n"
		 << "hx: " << hx << "\n"
		 << "hz: " << hz << "\n"
		 << "total time: " << ttotal << "\n"
		 << "dt: " << tstep << "\n"
		 << "meas_time: " << nmeas << "\n"
		 << "max_bond_dim: " << bonddim << endl;
		 
	return 0;
}
