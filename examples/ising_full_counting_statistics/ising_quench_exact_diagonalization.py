"""Exact diagonalization of the quench of ising_quench, to check the TN results (small N).

Usage: python3 ising_quench_exact_diagonalization.py input_ising_quench.txt
Output: data/<run>/entropy_exact_diagonalization.txt and gf_t<t>_exact_diagonalization.txt, compared with the TN files.
"""
import os
import sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from exact_diagonalization_tools import read_input, run_directory, save_and_compare, theta_grid
from ising_common_exact_diagonalization import ising_hamiltonian, generating_function_table


def entanglement_entropy(psi, b):
    """Von Neumann entropy (natural log) across the bond (b, b+1)."""
    s = np.linalg.svd(psi.reshape(2 ** b, -1), compute_uv=False) ** 2
    s = s[s > 1e-12]
    return -np.sum(s * np.log(s))


def main(input_file):
    p = read_input(input_file)
    N, J, hx, hz = p.get('N', 16), p.get('J', 1.), p.get('hx', 0.1), p.get('hz', 1.)
    T, max_dim, state = p.get('T', 5.), p.get('max_dim', 128), p.get('state', 'up')
    t_measure, number_points = p.get('t_measure', 0.5), p.get('number_points', 100)
    max_block_size = p.get('max_block_size', N // 2)

    # product state along x: '0' -> |+x>, '1' -> |-x>
    config = {'up': '0' * N, 'down': '1' * N, 'wall': '0' * (N // 2) + '1' * (N - N // 2)}[state]
    plus, minus = np.array([1., 1.]) / np.sqrt(2), np.array([1., -1.]) / np.sqrt(2)
    psi0 = np.array([1.])
    for c in config:
        psi0 = np.kron(psi0, plus if c == '0' else minus)

    E, V = np.linalg.eigh(ising_hamiltonian(N, J, hx, hz))
    theta = theta_grid(number_points)
    run_dir = run_directory('ising_quench_N%d_J%.2f_hx%.2f_hz%.2f_%s_T%g_dt%g_D%d' % (N, J, hx, hz, state, T, p.get('dt', 0.01), max_dim))
    gf_header = 'theta ' + ' '.join(f'ReG_{l} ImG_{l}' for l in range(1, max_block_size + 1))

    comparisons = {}
    entropy = []
    for t in np.arange(0, T + 1e-9, t_measure):
        psi = V @ (np.exp(-1j * E * t) * (V.conj().T @ psi0))
        entropy.append([t] + [entanglement_entropy(psi, b) for b in range(1, N)])
        comparisons['gf_t%.2f' % t] = save_and_compare(generating_function_table(psi, N, theta, max_block_size),
                                                       gf_header, '%sgf_t%.2f.txt' % (run_dir, t))
    comparisons['entropy'] = save_and_compare(entropy, 't ' + ' '.join(f'S_{b}' for b in range(1, N)), run_dir + 'entropy.txt')
    return comparisons


if __name__ == '__main__':
    main(sys.argv[1])
