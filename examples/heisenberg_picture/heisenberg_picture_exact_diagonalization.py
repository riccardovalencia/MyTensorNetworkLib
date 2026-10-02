"""Exact <psi_0| Z_j(t) |psi_0> for heisenberg_picture (dense, N <= ~12).

    H = sum_j (Jx X_j X_{j+1} + Jy Y_j Y_{j+1} + Jz Z_j Z_{j+1}) + sum_j (hx X_j + hy Y_j + hz Z_j)

Compares the exact result with the TN results in the Schroedinger and in the Heisenberg picture
("expectation"; the error is the one of the propagator, TEBD gates or first-order MPO) and reports
how much the two pictures differ ("pictures": they should agree up to truncation).

Usage: python3 heisenberg_picture_exact_diagonalization.py input.txt
"""
import os
import sys
from functools import reduce
import numpy as np
from scipy.linalg import expm

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from exact_diagonalization_tools import read_input, save_and_compare, tn_root  # noqa: E402

PAULI = {'x': np.array([[0., 1.], [1., 0.]], dtype=complex),
         'y': np.array([[0., -1j], [1j, 0.]]),
         'z': np.array([[1., 0.], [0., -1.]], dtype=complex)}
# single-site states of make_product_state (basis states 1 = up_z, 2 = down_z)
PRODUCT_STATES = {'z': {'0': [1, 0], '1': [0, 1]},
                  'x': {'0': [1, 1], '1': [1, -1]},
                  'y': {'0': [1, 1j], '1': [1, -1j]}}


def site_operator(o, j, N):
    """o on site j (1-indexed); site 1 is the leftmost tensor factor."""
    return reduce(np.kron, [o if k == j else np.eye(2) for k in range(1, N + 1)])


def hamiltonian(N, J, h):
    H = sum(J[d] * site_operator(PAULI[d], j, N) @ site_operator(PAULI[d], j + 1, N)
            for d in 'xyz' for j in range(1, N))
    return H + sum(h[d] * site_operator(PAULI[d], j, N) for d in 'xyz' for j in range(1, N + 1))


def product_state(config, basis):
    psi = reduce(np.kron, [np.array(PRODUCT_STATES[basis][c], dtype=complex) for c in config])
    return psi / np.linalg.norm(psi)


def main(input_file):
    p = read_input(input_file)
    N, site = p.get('N', 10), p.get('site', 5)
    J = {'x': p.get('Jx', 1.), 'y': p.get('Jy', 0.), 'z': p.get('Jz', 0.5)}
    h = {'x': p.get('hx', 0.4), 'y': p.get('hy', 0.3), 'z': p.get('hz', 0.6)}
    T, dt, t_measure = p.get('T', 2.), p.get('dt', 0.02), p.get('t_measure', 0.1)
    config, basis = str(p.get('config', '0110100101')), p.get('basis', 'x')

    H = hamiltonian(N, J, h)
    Z = site_operator(PAULI['z'], site, N)
    psi = product_state(config, basis)
    times = np.arange(0, T + 1e-9, t_measure)
    step = expm(-1j * t_measure * H)
    values = []
    for _ in times:
        values.append(np.vdot(psi, Z @ psi).real)
        psi = step @ psi

    root = tn_root('data/heisenberg_%s_N%d_site%d_J%.2f_%.2f_%.2f_h%.2f_%.2f_%.2f_dt%.4f', p.get('propagator', 'gates'),
                   N, site, J['x'], J['y'], J['z'], h['x'], h['y'], h['z'], dt)
    data = np.column_stack([times, values, values])
    comparisons = {'expectation': save_and_compare(root + '_ED.txt', data, 't Z_Schroedinger Z_Heisenberg', root + '.txt'),
                   'pictures': None}
    if os.path.exists(root + '.txt'):
        tn = np.loadtxt(root + '.txt', ndmin=2)
        comparisons['pictures'] = {'Z_S - Z_H': np.abs(tn[:, 1] - tn[:, 2]).max()}
        print(f'|Schroedinger - Heisenberg| = {comparisons["pictures"]["Z_S - Z_H"]:.2e}')
    return comparisons


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'input_heisenberg_picture.txt')
