#include <itensor/all.h>
#include <iostream>
#include <sstream> // for ostringstream
#include <vector>
#include <string>
#include <iomanip>
#include <complex>
#include "mytn.h"
#include <filesystem>

using namespace std;
using namespace itensor;
namespace fs = std::filesystem;	
	
//----------------------------------------------------------------------MAIN-----------------------------------------------------------------

int main(int argc , char* argv[])	{
	
	//input	
	if(argc != 9)
		{
		cout << "Input: " << argv[0] << " <N> <hxChoice> <hzChoice> <timestep> <nmeas> <numberPoints> <maxLength> <LOCALE (0) / CLUSTER (1)>" << endl;
		return -1;
		}
		
	int N, nmeas , numberPoints , maxLength , hxChoice , hzChoice , localvscluster;
	double tstep; 
	double wall_time = 12. * 60. * 60. ; //hours * minutes * seconds (wall time in seconds)
	
	// all files are read/written in the folder data/
	fs::create_directories("data");
	fs::current_path("data");

	get_data_meas( argv , &N , &hxChoice , &hzChoice , &tstep, &nmeas , &numberPoints , &maxLength , &localvscluster);
	
	//int numberStep = int( ttotal/tstep + (1E-9*(ttotal/tstep)) );
	
	cout << "Input:\n"
		 << "N : " << N << "\n"
		 << "tstep " << tstep << "\n"
		 << "nmeas : " << nmeas << "\n"
		 << "numberPoints : " << numberPoints << "\n"
		 << "maxLength : " << maxLength << "\n"
		 << "hxChoice : " << hxChoice << "\n"
		 << "hzChoice : " << hzChoice << endl;
	
	cout << "Getting the sites_file and psi_file" << endl;
	
	stringstream psi_file , sites_file , save_real , save_imag , saveRealMoments , saveImagMoments , out;

	build_file_full_counting( &sites_file , &psi_file , &save_real , &save_imag , &saveRealMoments , &saveImagMoments , N ); 
	 	
	cout << "The files found are called\n"
		 << "Sites -> " << sites_file.str() << "\n"
		 << "Psi -> " << psi_file.str() << "\n" << endl;
 	
	SpinHalf sites = readFromFile<SpinHalf>(sites_file.str());
	MPS psi;
		
	int n = 0;
	
	vector<double> Mx;
	
	cout << "Wall time : " << wall_time << endl;
	time_t c_begin = time(NULL);
	time_t c_step_start;
	time_t c_step_end;
	
	do
		{
		out.str("");
		out << psi_file.str() << n;
		if( fileExists( out.str() ) == true)
			{
			c_step_start = time(NULL);
			
			stringstream saveRealMomentsTwo;
			saveRealMomentsTwo << saveRealMoments.str() << n << ".dat";
			
			stringstream save_real_two;
			save_real_two << save_real.str() << n << ".dat";
			
			//if( fileExists( saveRealMomentsTwo.str() ) == false)
				//{
				//readFromFile(out.str(),psi);
				//measuring_moments( &psi , sites , N , n , maxLength , &saveRealMoments , &saveImagMoments);
							
				//if( n % 200 == 0 ) 
					//{
					//psi.position(N/2); 									//arbitrary, in order to avoid problems in normalize(psi)
					//normalize(psi);
					//measure_generating_function( &psi , sites , N , n  , maxLength , numberPoints , &save_real , &save_imag);
					//}
				//}
			//else if( fileExists( save_real_two.str() ) == false )
				//{
				//if( n % 200 == 0 ) 
					//{
					//readFromFile(out.str(),psi);
					//psi.position(N/2); 									//arbitrary, in order to avoid problems in normalize(psi)
					//normalize(psi);
					//measure_generating_function( &psi , sites , N , n  , maxLength , numberPoints , &save_real , &save_imag);
					//}		
				//}
			
			//if( n<= 200 )
				//{
				psi = readFromFile<MPS>(out.str(),sites);
				psi.position(N/2); 										//arbitrary, in order to avoid problems in normalize(psi)
				psi.normalize();
				measure_generating_function( &psi , sites , N , n  , maxLength , numberPoints , &save_real , &save_imag);
				//}
			
			c_step_end = time(NULL);
			time_t time_elapsed_step = c_step_end - c_step_start;
			time_t time_elapsed_total = c_step_end - c_begin;
			cout << "Single Step Time : " << time_elapsed_step << "\n"
				 << "Total Time : " <<  time_elapsed_total << endl;
			}
		n += nmeas;
		} while( fileExists( out.str() ) == true || n < 1000 );
	
	n -= 2*nmeas;
	out.str("");
	out << psi_file.str() << n;
	cout << "Last file opened is " << out.str() << endl;
	
	return 0;
}
