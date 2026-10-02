"""Exact diagonalization of the quench of ising_quench, to check the TN results (small N): evolved state,
thermal state with the energy of the initial state (beta fixed exactly; the TN program reaches it within
one imaginary-time step 2 dbeta) and the distance of the generating functions from the thermal ones.

Usage: python3 ising_quench_exact_diagonalization.py input_ising_quench.txt
Output (data/<run>/): <file>_exact_diagonalization.txt for xj.txt, entropy.txt, gf_t<t>.txt, thermal_gf.txt,
thermal_xj.txt and distance_to_thermal.txt, compared with the TN files.
"""
import os
import sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from exact_diagonalization_tools import read_input, run_directory, save_and_compare, theta_grid
from ising_common_exact_diagonalization import (X, ising_hamiltonian, generating_function_table, initial_state,
                                                site_operator, thermal_density_matrix)


def entanglement_entropy(psi, b):
    """Von Neumann entropy (natural log) across the bond (b, b+1)."""
    s = np.linalg.svd(psi.reshape(2 ** b, -1), compute_uv=False) ** 2
    s = s[s > 1e-12]
    return -np.sum(s * np.log(s))


def distances(table, reference):
    """max_theta |G_l - G_l^reference| for every block l, from two tables theta, Re G_1, Im G_1, ..."""
    G = table[:, 1::2] + 1j * table[:, 2::2]
    G_reference = reference[:, 1::2] + 1j * reference[:, 2::2]
    return np.abs(G - G_reference).max(axis=0)


def main(input_file):
    p = read_input(input_file)
    N, J, hx, hz = p.get('N', 16), p.get('J', 1.), p.get('hx', 0.1), p.get('hz', 1.)
    T, dt, max_dim, state = p.get('T', 5.), p.get('dt', 0.01), p.get('max_dim', 128), p.get('state', 'up')
    t_measure, number_points = p.get('t_measure', 0.5), p.get('number_points', 100)
    max_block_size, dbeta = p.get('max_block_size', N // 2), p.get('dbeta', 0.001)

    psi0 = initial_state(N, state)
    H = ising_hamiltonian(N, J, hx, hz)
    E, V = np.linalg.eigh(H)
    theta = theta_grid(number_points)
    run_dir = run_directory('ising_quench_N%d_J%.2f_hx%.2f_hz%.2f_%s_T%g_dt%g_D%d_dbeta%g' % (N, J, hx, hz, state, T, dt, max_dim, dbeta))
    gf_header = 'theta ' + ' '.join(f'ReG_{l} ImG_{l}' for l in range(1, max_block_size + 1))

    rho_thermal, beta = thermal_density_matrix(E, V, np.vdot(psi0, H @ psi0).real)
    print(f'thermal state: beta = {beta:.6f}')
    G_thermal = generating_function_table(rho_thermal, N, theta, max_block_size)
    Xj = [site_operator(X, j, N) for j in range(N)]
    comparisons = {'thermal_gf': save_and_compare(G_thermal, gf_header, run_dir + 'thermal_gf.txt'),
                   'thermal_xj': save_and_compare([[j + 1, np.trace(rho_thermal @ Xj[j]).real] for j in range(N)], 'j X_j',
                                                  run_dir + 'thermal_xj.txt')}

    xj, entropy, distance = [], [], []
    for t in np.arange(0, T + 1e-9, t_measure):
        psi = V @ (np.exp(-1j * E * t) * (V.conj().T @ psi0))
        xj.append([t] + [np.vdot(psi, x @ psi).real for x in Xj])
        entropy.append([t] + [entanglement_entropy(psi, b) for b in range(1, N)])
        G = generating_function_table(psi, N, theta, max_block_size)
        distance.append([t] + list(distances(G, G_thermal)))
        comparisons['gf_t%.2f' % t] = save_and_compare(G, gf_header, '%sgf_t%.2f.txt' % (run_dir, t))
    comparisons['xj'] = save_and_compare(xj, 't ' + ' '.join(f'X_{j}' for j in range(1, N + 1)), run_dir + 'xj.txt')
    comparisons['entropy'] = save_and_compare(entropy, 't ' + ' '.join(f'S_{b}' for b in range(1, N)), run_dir + 'entropy.txt')
    comparisons['distance'] = save_and_compare(distance, 't ' + ' '.join(f'D_{l}' for l in range(1, max_block_size + 1)),
                                               run_dir + 'distance_to_thermal.txt')
    return comparisons


if __name__ == '__main__':
    main(sys.argv[1])
