#ifndef H_SINGLE_PARTICLE_FERMIONS
#define H_SINGLE_PARTICLE_FERMIONS

#include <itensor/all.h>

using namespace std;
using namespace itensor;

ITensor
H_number_conserving_fermions(const int N , const vector<double> J, const vector<double> h, const bool spinful);

ITensor
H_number_conserving_fermions_impurity(const int N , const vector<double> J, const vector<double> h, const bool spinful);
#endif
