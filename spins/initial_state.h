#ifndef INITIAL_STATE
#define INITIAL_STATE

#include <itensor/all.h>

using namespace std;
using namespace itensor;

//----------------------------------------------------------------------
//all spins UP along X
void 
initial_state_all_UP( const SpinHalf , MPS* , const int );

//----------------------------------------------------------------------
//all spins DOWN along X
void
initial_state_all_DOWN( const SpinHalf , MPS* , const int );

//----------------------------------------------------------------------
//half chain DOWN along X, half chain UP along X (domain wall with single kink)
void
initial_state_DOMAIN_WALL( const SpinHalf , MPS* , const int );

#endif
