"""Exact ground state of the spin chain of spin_chain_ground_state (sparse Lanczos, N <= ~20).

    H = sum_j (Jx X_j X_{j+1} + Jy Y_j Y_{j+1} + Jz Z_j Z_{j+1})
      + sum_j (J2x X_j X_{j+2} + J2y Y_j Y_{j+2} + J2z Z_j Z_{j+2}) + sum_j (hx X_j + hz Z_j)

With conserve_sz = 1 the ground state is searched in the sector S^z = 0, as in the TN program.
Compares the energy (and reports the TN variance <H^2> - <H>^2, zero for an exact eigenstate) and
the profile j, <Z_j>, <Z_j Z_{j+1}>, entanglement entropy (log2) of the bond (j, j+1).

Usage: python3 spin_chain_ground_state_exact_diagonalization.py input.txt
"""
import os
import sys
from functools import reduce
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from exact_diagonalization_tools import read_input, run_directory, save_and_compare  # noqa: E402

PAULI = {'x': sp.csr_matrix([[0., 1.], [1., 0.]]),
         'y': sp.csr_matrix([[0., -1j], [1j, 0.]]),
         'z': sp.csr_matrix([[1., 0.], [0., -1.]])}


def site_operator(o, j, N):
    """Sparse operator o on site j (0-indexed); site 0 is the leftmost tensor factor."""
    return reduce(lambda a, b: sp.kron(a, b, format='csr'),
                  [o if k == j else sp.identity(2, format='csr') for k in range(N)])


def spin_chain_hamiltonian(N, J, J2, h):
    """J, J2, h: dicts direction -> coupling (nearest, next-nearest neighbours, fields)."""
    paulis = {d: [site_operator(PAULI[d], j, N) for j in range(N)] for d in 'xyz'}
    H = sp.csr_matrix((2 ** N, 2 ** N), dtype=complex)
    for d in 'xyz':
        for distance, couplings in ((1, J), (2, J2)):
            if couplings[d] != 0:
                H = H + couplings[d] * sum(paulis[d][j] @ paulis[d][j + distance] for j in range(N - distance))
        if h.get(d, 0) != 0:
            H = H + h[d] * sum(paulis[d])
    return H


def zero_magnetization_sector(N):
    """Indices of the basis states (bit 1 = spin down) with as many up as down spins."""
    return np.array([s for s in range(2 ** N) if bin(s).count('1') == N // 2])


def ground_state(H, N, conserve_sz):
    """Lowest eigenvalue and normalized eigenvector of H (in the S^z = 0 sector if conserve_sz)."""
    if not conserve_sz:
        energies, vectors = eigsh(H, k=1, which='SA', tol=1e-12)
        return energies[0], vectors[:, 0]
    sector = zero_magnetization_sector(N)
    energies, vectors = eigsh(H[sector][:, sector], k=1, which='SA', tol=1e-12)
    psi = np.zeros(2 ** N, dtype=complex)
    psi[sector] = vectors[:, 0]
    return energies[0], psi


def entanglement_entropy(psi, N, j):
    """Von Neumann entropy (log2) of sites 1..j (1-indexed)."""
    s = np.linalg.svd(psi.reshape(2 ** j, 2 ** (N - j)), compute_uv=False) ** 2
    s = s[s > 1e-16]
    return float(-np.sum(s * np.log2(s)))


def profile(psi, N):
    """Rows j, <Z_j>, <Z_j Z_{j+1}>, S(j, j+1) (last two zero for j = N)."""
    z = [site_operator(PAULI['z'], j, N) for j in range(N)]
    rows = []
    for j in range(N):
        zz = np.vdot(psi, z[j] @ (z[j + 1] @ psi)).real if j < N - 1 else 0.
        entropy = entanglement_entropy(psi, N, j + 1) if j < N - 1 else 0.
        rows.append([j + 1, np.vdot(psi, z[j] @ psi).real, zz, entropy])
    return np.array(rows)


def main(input_file):
    p = read_input(input_file)
    N = p.get('N', 16)
    J = {d: p.get(f'J{d}', 1.) for d in 'xyz'}
    J2 = {d: p.get(f'J2{d}', 0.) for d in 'xyz'}
    h = {'x': p.get('hx', 0.), 'z': p.get('hz', 0.)}
    conserve_sz = p.get('conserve_sz', 0) == 1

    energy, psi = ground_state(spin_chain_hamiltonian(N, J, J2, h), N, conserve_sz)

    run = 'spin_chain_N%d_J%.3f_%.3f_%.3f_J2%.3f_%.3f_%.3f_hx%.3f_hz%.3f' % (
        N, J['x'], J['y'], J['z'], J2['x'], J2['y'], J2['z'], h['x'], h['z'])
    run_dir = run_directory(run + '_sz' if conserve_sz else run)
    comparisons = {'energy': None}
    print(f'ED ground-state energy: {energy:.12f}')
    if os.path.exists(run_dir + 'energy.txt'):
        tn_energy, tn_variance, _ = np.loadtxt(run_dir + 'energy.txt')
        comparisons['energy'] = {'E': abs(tn_energy - energy), 'variance': abs(tn_variance)}
        print(f'Comparison with {run_dir}energy.txt: |E_TN - E_ED| = {abs(tn_energy - energy):.2e}, TN variance = {tn_variance:.2e}')
    comparisons['profile'] = save_and_compare(profile(psi, N), 'j <Z_j> <Z_jZ_j+1> S', run_dir + 'profile.txt')
    return comparisons


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'input_spin_chain_ground_state.txt')
