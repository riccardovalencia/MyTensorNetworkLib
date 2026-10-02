"""Exact thermal state of ising_thermal, to check the TN results (small N).

beta is fixed by requiring the energy density of |+x...+x>, -J ((N-1)/N + hx); the TN program reaches
it with steps of 2 dbeta, so small differences of order dbeta are expected.

Usage: python3 ising_thermal_exact_diagonalization.py input_ising_thermal.txt
Output: data/<run>/gf_exact_diagonalization.txt, compared with the TN file gf.txt.
"""
import os
import sys
import numpy as np
from scipy.optimize import brentq
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from exact_diagonalization_tools import read_input, run_directory, save_and_compare, theta_grid
from ising_common_exact_diagonalization import ising_hamiltonian, generating_function_table


def main(input_file):
    p = read_input(input_file)
    N, J, hx, hz = p.get('N', 16), p.get('J', 1.), p.get('hx', 0.1), p.get('hz', 1.)
    number_points, max_block_size = p.get('number_points', 100), p.get('max_block_size', N // 2)

    E, V = np.linalg.eigh(ising_hamiltonian(N, J, hx, hz))
    energy_target = -J * ((N - 1.) / N + hx)

    def energy_density(beta):
        w = np.exp(-beta * (E - E.min()))
        return np.sum(w * E) / np.sum(w) / N

    beta = brentq(lambda b: energy_density(b) - energy_target, 0., 50.)
    w = np.exp(-beta * (E - E.min()))
    rho = V @ np.diag(w / w.sum()) @ V.conj().T
    print(f'beta = {beta:.6f}')

    run_dir = run_directory('ising_thermal_N%d_J%.2f_hx%.2f_hz%.2f_dbeta%g_D%d' % (N, J, hx, hz, p.get('dbeta', 0.001), p.get('max_dim', 1000)))
    return {'gf': save_and_compare(generating_function_table(rho, N, theta_grid(number_points), max_block_size),
                                   'theta ' + ' '.join(f'ReG_{l} ImG_{l}' for l in range(1, max_block_size + 1)), run_dir + 'gf.txt')}


if __name__ == '__main__':
    main(sys.argv[1])
