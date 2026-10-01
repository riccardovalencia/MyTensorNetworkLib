#ifndef TEBD_TESTING_H
#define TEBD_TESTING_H

#include <itensor/all.h>
#include "MyClasses.h"

using namespace std;
using namespace itensor;

vector<MyBondGate>
gates_rydberg_up_to_VNN_testing(const SiteSet sites , const vector<double> Deltaj, const vector<double> Omegaj, const vector<double> Vj, const double dt);


vector<MyBondGate>
gates_rydberg_up_to_VNNN_testing(const SiteSet sites , const vector<double> Deltaj, const vector<double> Omegaj, const vector<double> Vj, const double dt);

vector<MyBondGate>
gates_rydberg_up_to_VNNN_testing_2(const SiteSet sites , const vector<double> Deltaj, const vector<double> Omegaj, const vector<double> Vj, const double dt);

vector<MyBondGate>
gates_spin_local_field_testing(const SiteSet sites , vector<double> omegaj, const double dt);


vector<MyBondGateDiss>
gates_local_lindbland_testing(const SiteSet sites , vector<ITensor> Lj, vector<double> gammaj , const double dt);

vector<MyBondGateDiss>
gates_nearest_neighbour_local_lindbland_testing(const SiteSet sites , vector<MyTrainITensor> TTrain, const double dt);

#endif
