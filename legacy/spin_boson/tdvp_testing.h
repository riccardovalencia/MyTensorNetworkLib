#ifndef TEBD_LONG_RANGE_H
#define TEBD_LONG_RANGE_H

#include <itensor/all.h>
#include "MyClasses.h"
using namespace std;
using namespace itensor;


MPS
TDVP_lindbland_time_evolve(MPS psi_t, MPO H , vector<MyBondGateDiss> gates_D , Args TEBD_args, bool dissipative , double dt , double T , int steps_save_state, bool normalize, string file_root, double t_start = 0.);

MPS
TDVP_split_lindbland_time_evolve(MPS psi_t, MPO Hbra , MPO Hket, vector<MyBondGateDiss> gates_D , Args TEBD_args, bool dissipative , double dt , double T , int steps_save_state, bool normalize, string file_root, double t_start = 0.);

MPS
TDVP_time_evolve(MPS psi_t,SiteSet sites, MPO H , Args TDVP_args, double dt , double T , int steps_save_state, bool normalize, string file_root, double t_start = 0.);

#endif
