#ifndef MYTN_SPINS_GET_DATA_H
#define MYTN_SPINS_GET_DATA_H

// Command-line parsing for the Ising-chain drivers of the spins module.

using namespace std;

// argv = state N J hxChoice hzChoice ttotal tstep nmeas bonddim localvscluster;
// hx and hz are picked from fixed lists by their index.
void
get_data( char* argv[] , int *state , int *N , double *J , double *hx , double *hz , double *ttotal , double *tstep , int *nmeas , int *bonddim , int *localvscluster );

// argv = N hxChoice hzChoice tstep nmeas numberPoints maxLength localvscluster.
void
get_data_meas( char* argv[] , int *N , int *hxChoice , int *hzChoice , double *tstep , int *nmeas , int *numberPoints , int *maxLength , int *localvscluster );

// argv = N tstep nmeas localvscluster.
void
get_data_entropy( char* argv[] , int *N , double *tstep , int *nmeas , int *localvscluster );

#endif
