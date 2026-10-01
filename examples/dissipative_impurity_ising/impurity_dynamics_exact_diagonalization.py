"""Exact solution of the integrable case of impurity_dynamics (Jxx = Jzzz = 0), to check the TN results.

    H = -sum_j Z_j Z_{j+1} + hx sum_j X_j,   L = sqrt(gamma) Z_1.
After a Jordan-Wigner transformation (X_j = 1 - 2 c^dag_j c_j) H is quadratic and Z_1 is a single
Majorana operator, so the Lindblad dynamics closes on the 2N x 2N Majorana covariance matrix C:
this gives the exact dynamics for large N (no exact diagonalization of the many-body Hamiltonian).
Starting from the ground state of H, the script computes
  - the profile <X_j>(t) up to Tness,
  - the autocorrelation <Z_1(t) Z_1(0)> in the state reached at Tness, up to time T.

Usage: python3 impurity_dynamics_exact_diagonalization.py input_impurity_dynamics.txt
Output: data/<root>_exact_diagonalization_xj.txt and ..._z1z1.txt, compared with the TN files.
"""
import os
import sys
import numpy as np
from scipy.integrate import solve_ivp
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from exact_diagonalization_tools import read_input, save_and_compare

p = read_input(sys.argv[1])
N, hx, gamma = p.get('N', 10), p.get('hx', 0.5), p.get('gamma', 0.5)
Jxx, Jzzz = p.get('Jxx', 0.2), p.get('Jzzz', 0.)
Tness, T, dt, max_dim = p.get('Tness', 5.), p.get('T', 5.), p.get('dt', 0.05), p.get('max_dim', 128)
if Jxx != 0 or Jzzz != 0:
    sys.exit('The free-fermion solution requires Jxx = Jzzz = 0 (integrable case).')
J = 1.
t_measure, t_corr = p.get('t_measure', 0.2), p.get('t_corr', 0.05)   # as in the TN program


def bdg_hamiltonian(N, J, h):
    """Bogoliubov-de Gennes matrix of the Ising chain in the Dirac-fermion basis (c, c^dag)."""
    A = np.diag(np.full(N, h)) + np.diag(np.full(N - 1, -J / 2), 1) + np.diag(np.full(N - 1, -J / 2), -1)
    B = np.diag(np.full(N - 1, -J / 2), 1) + np.diag(np.full(N - 1, J / 2), -1)
    return np.block([[A, B], [-B, -A]])


def majorana_hamiltonian(N, J, h):
    """Single-particle Hamiltonian of the Ising chain in the Majorana basis."""
    block = np.diag(np.full(N, -h / 2)) + np.diag(np.full(N - 1, J / 2), -1)
    H = np.zeros((2 * N, 2 * N), dtype=complex)
    H[:N, N:] = block
    H[N:, :N] = -block.T
    return 1j * H


def ground_state_covariance(N, J, h):
    """Majorana covariance matrix of the BdG vacuum (ground state)."""
    _, U = np.linalg.eigh(bdg_hamiltonian(N, J, h))
    U = U[:, ::-1]
    C_dirac = U @ np.diag(np.r_[np.ones(N), np.zeros(N)]) @ U.conj().T
    Om = np.block([[np.eye(N), np.eye(N)], [1j * np.eye(N), -1j * np.eye(N)]])
    return Om @ C_dirac @ Om.conj().T


H = majorana_hamiltonian(N, J, hx)
P1 = np.zeros((2 * N, 2 * N))
P1[0, 0] = 1.      # dephasing acts on the first Majorana operator (Z_1)


def covariance_rhs(t, c):
    """Lindblad dynamics of the covariance matrix."""
    C = c.reshape(H.shape)
    dC = -4j * (H @ C - C @ H) - 2 * gamma * (P1 @ C + C @ P1) + 4 * gamma * P1 @ C @ P1
    return dC.ravel()


def autocorrelation_rhs(t, c):
    """Regression-theorem dynamics for <Z_1(t) Z_1(0)>: the dephasing damps every Majorana
    operator except the first one at rate 2 gamma."""
    C = c.reshape(H.shape)
    damping = np.eye(H.shape[0]) - P1
    return ((-4j * H - 2 * gamma * damping) @ C).ravel()


def integrate(rhs, C0, times):
    sol = solve_ivp(rhs, (0, times[-1]), C0.ravel().astype(complex), t_eval=times, method='DOP853', rtol=1e-12, atol=1e-12)
    return [sol.y[:, k].reshape(C0.shape) for k in range(len(times))]


# relaxation: <X_j> = Im C[j, N+j]
times = np.arange(t_measure, Tness + 1e-9, t_measure)
covariances = integrate(covariance_rhs, ground_state_covariance(N, J, hx), np.r_[times, Tness] if times[-1] < Tness - 1e-9 else times)
xj = [[t] + list(np.imag(np.diag(C[:N, N:]))) for t, C in zip(times, covariances)]

# autocorrelation from the covariance matrix at Tness
times_corr = np.arange(t_corr, T + 1e-9, t_corr)
z1z1 = [[t, C[0, 0].real, C[0, 0].imag, abs(C[0, 0])] for t, C in zip(times_corr, integrate(autocorrelation_rhs, covariances[-1], times_corr))]

root = 'data/impurity_N%d_Jxx%.3f_Jzzz%.3f_hx%.3f_gamma%.3f_dt%.4f_D%d' % (N, Jxx, Jzzz, hx, gamma, dt, max_dim)
save_and_compare(root + '_exact_diagonalization_xj.txt', xj, 't ' + ' '.join(f'X_{j}' for j in range(1, N + 1)), root + '_xj.txt')
save_and_compare(root + '_exact_diagonalization_z1z1.txt', z1z1, 't Re Im abs', root + '_Tness%.1f_z1z1.txt' % Tness)
