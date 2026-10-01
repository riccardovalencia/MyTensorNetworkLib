/**
 * @file input.h
 * @brief Data management: command-line parsing for the drivers (argv[0] is the program name).
 */
#ifndef MYTN_IO_INPUT_H
#define MYTN_IO_INPUT_H

using namespace std;

/** @name Ising chain (spins) */
///@{

/**
 * @brief argv = state N J hxChoice hzChoice ttotal tstep nmeas bonddim localvscluster.
 *
 * hx and hz are taken from fixed lists by index:
 * hx in {0, 0.1, 0.2, 0.4, -0.1, -0.2, -0.4}, hz in {0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 3}.
 */
void
get_data( char* argv[] , int *state , int *N , double *J , double *hx , double *hz , double *ttotal , double *tstep , int *nmeas , int *bonddim , int *localvscluster );

/** @brief argv = N hxChoice hzChoice tstep nmeas numberPoints maxLength localvscluster. */
void
get_data_meas( char* argv[] , int *N , int *hxChoice , int *hzChoice , double *tstep , int *nmeas , int *numberPoints , int *maxLength , int *localvscluster );

/** @brief argv = N tstep nmeas localvscluster. */
void
get_data_entropy( char* argv[] , int *N , double *tstep , int *nmeas , int *localvscluster );
///@}

/** @name Bosonic quantum east model */
///@{

/** @brief argv[1..6] = size s c n0 symmetry_sector cut_off_fock_space. */
void
get_data_system( char* argv[] , int *size , double *s , double *c, int *n0 , double *simmetry_sector , int *cut_off_fock_space );

/** @brief argv[7..11] = bond_dimension lower_bound_singular_values scaling_bond_dimension precision_dmrg max_bond_dimension. */
void
get_data_DMRG( char* argv[] , int *bond_dimension, double *lower_bound_singular_values, double *scaling_bond_dimension, double *precision_dmrg, int *max_bond_dimension);

/**
 * @brief argv[6..10] = max_bond_dimension lower_bound_singular_values total_time delta_t number_steps_samplig.
 */
void
get_data_TEBD( char* argv[] , int *max_bond_dimension, double *lower_bound_singular_values, double *total_time, double *delta_t, int *number_steps_samplig);
///@}

#endif
