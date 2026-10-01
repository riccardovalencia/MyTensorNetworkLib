/**
 * @file adiabatic.h
 * @brief Adiabatic preparation of bosonic quantum east model (bosonic east model) states: the facilitated hopping
 *        J is ramped from 0 to e^{-s} with TEBD (see models/bosonic_east_model.h).
 */
#ifndef MYTN_DYNAMICS_ADIABATIC_H
#define MYTN_DYNAMICS_ADIABATIC_H

#include <itensor/all.h>

using namespace std;
using namespace itensor;

/**
 * @brief Linear ramp J(t) = e^{-s} t / T, evolving *psi_start in place until J reaches e^{-s}.
 * @param psi_start State to evolve.
 * @param sites     Boson site set.
 * @param s, c      Target bosonic east model parameters.
 * @param dt        Time step.
 * @param T         Duration of the ramp.
 * @param args      Truncation of the TEBD steps ("Cutoff", "MaxDim").
 */
void
evolve_adiabatic_linear_ramp( MPS *psi_start, const SiteSet sites, const double s, const double c, const double dt, const double T, const Args& args = Args("Cutoff=",1E-10,"MaxDim=",50));

/**
 * @brief Ramp J(t) = e^{-s} tanh(t / T), evolving until J is within a relative 1e-6 of e^{-s}.
 * @param T    Time scale of the ramp.
 * @param args Truncation of the TEBD steps ("Cutoff", "MaxDim").
 */
void
evolve_adiabatic_tanh_ramp( MPS *psi_start, const SiteSet sites, const double s, const double c, const double dt, const double T, const Args& args = Args("Cutoff=",1E-16,"MaxDim=",1000));

#endif
