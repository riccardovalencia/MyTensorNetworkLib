#ifndef BUILD_HAMILTONIAN_H
#define BUILD_HAMILTONIAN_H

#include <itensor/all.h>

using namespace itensor;
using namespace std;

// Hamiltonian H = - 0.5 n_0 ( exp(-s) \sigma_1^x -1) - 0.5 \sum_{j=1}^{L-1} n_j (exp(-s) \sigma_{j+1}^x - (1-2c) n_{j+1} -1) + 0.5 n_L * sector + 0.5 * n_L

MPO
H_mmGcbQEM( const SiteSet sites, int size , int n0, double symmetry , double s, double c);


MPO
H_mmGcbQEM_with_drift( const SiteSet sites, int size , int n0, double symmetry , double s, double c);

// Hamiltonian -H = 0.5 n_0 ( exp(-s) \sigma_1^x -1) + 0.5 \sum_{j=1}^{L-1} n_j (exp(-s) \sigma_{j+1}^x - (1-2c) n_{j+1} -1) + 0.5 n_L * sector 

MPO
H_mmGcbQEM_minus( const SiteSet sites, int size , int n0, double symmetry , double s, double c);

// Hamiltonian H = - 0.5 \sum_{j=1}^{L-1} n_j (exp(-s) \sigma_{j+1}^x - (1-2c) n_{j+1} -1) + \epsilon/2 \sum_{j=1}^L n_j n_j - t \sum_{j=1}^{L-1} (a_j^\dagger a_{j+1}  + h.c. )

MPO
H_mmGcbQEM_onsite_hopping( const SiteSet sites, int size , double s, double c , double epsilon, double t);


// Hamiltonian H = - 0.5 n_0 ( exp(-s) \sigma_1^x -1) - 0.5 \sum_{j=1}^{L-1} n_j (exp(-s) \sigma_{j+1}^x - (1-2c) n_{j+1} -1) + \epsilon/2 \sum_{j=1}^L n_j n_j - 0.5 n_L * sector 

MPO
H_mmGcbQEM_onsite( const SiteSet sites, int size , int n0, double symmetry , double s, double c , double epsilon);

// Hamiltonian H = - 0.5 n_0 ( exp(-s) \sigma_1^x -1) - 0.5 \sum_{j=1}^{L-1} n_j (exp(-s) \sigma_{j+1}^x -1) + c/2 \sum_{j=1}^L n_j n_j - 0.5 n_L * (exp(-s) * sector - 1) 

MPO
H_mmGcbQEM_onsite_nonext( const SiteSet sites, int size , int n0, double symmetry , double s, double c );


// ------------------------------------------------------------
// Build Hamiltonian of Bosonic quantum east model with next-neighbour density-density interaction.
// No on-site density-density interaction
// site 0 is not touched. This is used to find the ground state within a certain symmetry sector fixed by
// n0 from a site of the form |n0>|psi>.

MPO
H_mmGcbQEM_n0_untouched( const SiteSet sites, int size , int n0, double symmetry , double s, double c);


// ------------------------------------------------------------
// Build Hamiltonian of Bosonic quantum east model with next-neighbour density-density interaction.
// No on-site density-density interaction
// Notice that n0 is not a number, but its defined in the Hamiltonian itself. So, it is not a conserved quantity
// and we do not expect DMRG to conserve it (it will try to minimize the Hamiltonian and that's it. At most it 
// will be stuck for a long time if you initialize DMRG with a state |n_0> \otimes |\psi> since |n_0> is an eigenstate
// of the operators acting on the 0-th site). 

MPO
H_mmGcbQEM_n0_not_fixed( const SiteSet sites, int size , double symmetry , double s, double c);


// Exponential of the Hamiltonian. This is done in order to apply it to operators

MPO
exp_H_mmGcbQEM_n0_notfixed( const SiteSet sites, int size , int n0, double symmetry , double J, double c, double dt);


#endif
