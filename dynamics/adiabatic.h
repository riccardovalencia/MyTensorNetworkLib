/**
 * @file adiabatic.h
 * @brief Adiabatic preparation of bosonic quantum east model (bosonic east model) states: the facilitated hopping
 *        J is ramped from 0 to J_target with TEBD (see models/bosonic_east_model.h).
 */
#ifndef MYTN_DYNAMICS_ADIABATIC_H
#define MYTN_DYNAMICS_ADIABATIC_H

#include <itensor/all.h>

using namespace std;
using namespace itensor;

/**
 * @brief Linear ramp J(t) = J_target t / T, evolving *psi_start in place until J reaches J_target.
 * @param psi_start State to evolve.
 * @param sites     Boson site set.
 * @param J_target  Final facilitated hopping amplitude (e^{-s}).
 * @param U         Density-density coefficient (1 - 2c), constant during the ramp.
 * @param dt        Time step.
 * @param T         Duration of the ramp.
 * @param args      Truncation of the TEBD steps ("Cutoff", "MaxDim").
 */
void
evolve_adiabatic_linear_ramp( MPS *psi_start, const SiteSet sites, const double J_target, const double U, const double dt, const double T, const Args& args = Args("Cutoff=",1E-10,"MaxDim=",50));

/**
 * @brief Ramp J(t) = J_target tanh(t / T), evolving *psi_start in place until J is within a
 *        relative 1e-6 of J_target.
 * @param psi_start State to evolve.
 * @param sites     Boson site set.
 * @param J_target  Final facilitated hopping amplitude (e^{-s}).
 * @param U         Density-density coefficient (1 - 2c), constant during the ramp.
 * @param dt        Time step.
 * @param T         Time scale of the ramp.
 * @param args      Truncation of the TEBD steps ("Cutoff", "MaxDim").
 */
void
evolve_adiabatic_tanh_ramp( MPS *psi_start, const SiteSet sites, const double J_target, const double U, const double dt, const double T, const Args& args = Args("Cutoff=",1E-16,"MaxDim=",1000));

#endif
