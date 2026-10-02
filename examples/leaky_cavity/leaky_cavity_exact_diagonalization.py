"""Exact diagonalization of the leaky-cavity model of leaky_cavity, to check the TN results.

    H = omega0 a^dag a + h sum_j Z_j + V sum_j n_j n_{j+1} + (g/sqrt(N)) * coupling,
    coupling = (a + a^dag) sum_j X_j ("dicke") or sum_j (sigma^+_j a + sigma^-_j a^dag) ("tavis"),
    d rho/dt = -i[H, rho] + kappa D[a] rho,

with g = g_ratio * g_c, g_c = sqrt((|h| - V)(omega0^2 + kappa^2/4)/(2 omega0)), the boson truncated
at max_occ, initial state |0> (x) (cos(theta/2)|up_z> + sin(theta/2)|down_z>)^N.
The Lindblad equation is integrated for the full density matrix: only small N and max_occ.

Usage: python3 leaky_cavity_exact_diagonalization.py input_leaky_cavity.txt
Output: data/<run>/observables_exact_diagonalization.txt (t, <X_1>, <Z_1>, <a^dag a>),
        compared with the TN file observables.txt.
"""
import os
import sys
import numpy as np
import quimb as qu
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from exact_diagonalization_tools import read_input, run_directory, save_and_compare, lindblad_evolution


def main(input_file):
    p = read_input(input_file)
    N, max_occ = p.get('N', 3), p.get('max_occ', 2)
    h, g_ratio, V, kappa = p.get('h', 1.), p.get('g', 1.6), p.get('V', 0.), p.get('kappa', 1.)
    T, max_dim, coupling = p.get('T', 15.), p.get('max_dim', 1024), p.get('coupling', 'dicke')
    omega0, theta = 1., p.get('theta', 0.9) * np.pi
    t_measure = 0.1     # ED output interval (the rows matching the TN times are compared)

    g = g_ratio * np.sqrt(0.5 * (abs(h) - V) * (omega0**2 + kappa**2 / 4) / omega0)

    # operators: site 0 = boson, sites 1..N = spins (|up_z> = basis state 0)
    dims = [max_occ + 1] + [2] * N
    a = qu.qu(np.diag(np.sqrt(np.arange(1, max_occ + 1)), k=1), qtype='dop')
    X, Z = qu.pauli('X'), qu.pauli('Z')
    n = qu.qu(np.diag([0., 1.]), qtype='dop')
    sigma_minus = qu.qu([[0, 0], [1, 0]], qtype='dop')   # |up_z> -> |down_z>

    def op(o, j):
        return qu.ikron(o, dims, j)

    H = omega0 * op(a.H @ a, 0)
    for j in range(1, N + 1):
        H = H + h * op(Z, j)
    for j in range(1, N):
        H = H + V * qu.ikron([n, n], dims, [j, j + 1])
    for j in range(1, N + 1):
        if coupling == 'dicke':
            H = H + g / np.sqrt(N) * qu.ikron([a + a.H, X], dims, [0, j])
        else:
            H = H + g / np.sqrt(N) * (qu.ikron([a, sigma_minus.H], dims, [0, j]) + qu.ikron([a.H, sigma_minus], dims, [0, j]))

    spin = np.cos(theta / 2) * qu.basis_vec(0, 2) + np.sin(theta / 2) * qu.basis_vec(1, 2)
    psi0 = qu.kron(qu.basis_vec(0, max_occ + 1), *[spin] * N)
    rho0 = psi0 @ psi0.H

    times = np.arange(t_measure, T + 1e-9, t_measure)
    rhos = lindblad_evolution(rho0, H, [np.sqrt(kappa) * op(a, 0)], times)
    obs = [op(X, 1), op(Z, 1), op(a.H @ a, 0)]
    data = [[t] + [np.real(np.trace(O @ rho)) / np.real(np.trace(rho)) for O in obs] for t, rho in zip(times, rhos)]

    run = 'leaky_cavity_%s_N%d_maxocc%d_h%.2f_gratio%.2f_V%.2f_kappa%.2f_theta%g_T%g_dt%g_D%d' % (
        coupling, N, max_occ, h, g_ratio, V, kappa, p.get('theta', 0.9), T, p.get('dt', 0.01), max_dim)
    return {'obs': save_and_compare(data, 't X_1 Z_1 n_photon', run_directory(run) + 'observables.txt', [2, 3, 4])}


if __name__ == '__main__':
    main(sys.argv[1])
