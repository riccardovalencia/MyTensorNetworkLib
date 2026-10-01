"""Exact thermal state of ising_thermal, to check the TN results (small N).

beta is fixed by requiring the energy density of |+x...+x>, -J ((N-1)/N + hx); the TN program reaches
it with steps of 2 dbeta, so small differences of order dbeta are expected.

Usage: python3 ising_thermal_exact_diagonalization.py input_ising_thermal.txt
Output: data/<root>_exact_diagonalization_gf.txt, compared with the TN file <root>_gf.txt.
"""
import os
import sys
import numpy as np
from scipy.optimize import brentq
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from exact_diagonalization_tools import read_input, save_and_compare, theta_grid
from ising_common_exact_diagonalization import ising_hamiltonian, generating_function_table

p = read_input(sys.argv[1])
N, J, hx, hz = p.get('N', 16), p.get('J', 1.), p.get('hx', 0.1), p.get('hz', 1.)
number_points = 100     # as in the TN program

E, V = np.linalg.eigh(ising_hamiltonian(N, J, hx, hz))
energy_target = -J * ((N - 1.) / N + hx)


def energy_density(beta):
    w = np.exp(-beta * (E - E.min()))
    return np.sum(w * E) / np.sum(w) / N


beta = brentq(lambda b: energy_density(b) - energy_target, 0., 50.)
w = np.exp(-beta * (E - E.min()))
rho = V @ np.diag(w / w.sum()) @ V.conj().T
print(f'beta = {beta:.6f}')

root = 'data/ising_thermal_N%d_J%.2f_hx%.2f_hz%.2f' % (N, J, hx, hz)
save_and_compare(root + '_exact_diagonalization_gf.txt', generating_function_table(rho, N, theta_grid(number_points)),
                 'theta ' + ' '.join(f'ReG_{l} ImG_{l}' for l in range(1, N // 2 + 1)), root + '_gf.txt')
