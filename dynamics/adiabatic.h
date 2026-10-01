/**
 * @file adiabatic.h
 * @brief Adiabatic preparation of bosonic quantum east model (bQEM) states: the facilitated hopping
 *        J is ramped from 0 to e^{-s} with TEBD (see models/bqem.h).
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
 * @param s, c      Target bQEM parameters.
 * @param dt        Time step.
 * @param T         Duration of the ramp.
 */
void
adiabatic_transformation_linear_protocol( MPS *psi_start, const SiteSet sites, const double s, const double c, const double dt, const double T);

/**
 * @brief Ramp J(t) = e^{-s} tanh(t / T), evolving until J is within a relative 1e-6 of e^{-s}.
 * @param T Time scale of the ramp.
 */
void
adiabatic_transformation_tanh_protocol( MPS *psi_start, const SiteSet sites, const double s, const double c, const double dt, const double T);

#endif
