#ifndef TEBD_LONG_RANGE_H
#define TEBD_LONG_RANGE_H

#include <itensor/all.h>
#include "MyClasses.h"
using namespace std;
using namespace itensor;

vector<BondGate>
gates_tavis_cummings(const SiteSet sites , const double omega0 , const double h , const double g,const double dt,string photon_or_matter);

vector<BondGate>
gates_photon_matter(const SiteSet sites , const double omega0 , const double h , const double g, const double dt, string matter_or_photon, string type_of_coupling = "tavis", const double V = 0, string interaction_axis="z");

MPS
MPO_lindbland_time_evolve(MPS psi_t, MPO H , vector<MyBondGateDiss> gates_D , Args TEBD_args, bool dissipative , double dt , double T , int steps_save_state, bool normalize, string file_root, double t_start = 0.);

MPS
TEBD_long_range_int_lindbland_time_evolve(MPS psi_t, vector<BondGate> gates_H, vector<MyBondGateDiss> gates_D , Args TEBD_args, bool dissipative , double dt , double T , int steps_save_state, bool normalize, string file_root, double t_start = 0.);

#endif
