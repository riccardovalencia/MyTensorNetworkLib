"""Exact solution of impurity_dynamics, to check the TN results.

    H = -sum_j Z_j Z_{j+1} + Jxx sum_j X_j X_{j+1} + Jzzz sum_j Z_j Z_{j+2} + hx sum_j X_j,
    L = sqrt(gamma) Z_1.

Starting from the ground state of H, the script computes
  - the profile <X_j>(t) up to Tness,
  - the autocorrelation <Z_1(t) Z_1(0)> in the state reached at Tness, up to time T
    (quantum regression theorem: Z_1 rho is evolved with the Lindblad equation).

Integrable case (Jxx = Jzzz = 0): after a Jordan-Wigner transformation (X_j = 1 - 2 c^dag_j c_j) H is
quadratic and Z_1 is a single Majorana operator, so the Lindblad dynamics closes on the 2N x 2N
Majorana covariance matrix C: this gives the exact dynamics for large N.
Otherwise the Lindblad equation is integrated for the full density matrix (N <= 8).

Usage: python3 impurity_dynamics_exact_diagonalization.py input_impurity_dynamics.txt
Output: data/<run>/xj_exact_diagonalization.txt and z1z1_exact_diagonalization.txt, compared with the TN files.
"""
import os
import sys
from functools import reduce
import numpy as np
from scipy.integrate import solve_ivp
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from exact_diagonalization_tools import read_input, run_directory, save_and_compare, lindblad_evolution

MAX_DENSE_SITES = 8


# ----------------------------------------------------------
# integrable case: Majorana covariance matrix

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


def integrate(rhs, C0, times):
    sol = solve_ivp(rhs, (0, times[-1]), C0.ravel().astype(complex), t_eval=times, method='DOP853', rtol=1e-12, atol=1e-12)
    return [sol.y[:, k].reshape(C0.shape) for k in range(len(times))]


def free_fermion_solution(N, J, hx, gamma, times, times_corr):
    """Rows (t, <X_1>, ..., <X_N>) at times (the last one is Tness) and (t, Re, Im, abs) of the
    autocorrelation at times_corr."""
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

    covariances = integrate(covariance_rhs, ground_state_covariance(N, J, hx), times)
    xj = [[t] + list(np.imag(np.diag(C[:N, N:]))) for t, C in zip(times, covariances)]
    z1z1 = [[t, C[0, 0].real, C[0, 0].imag, abs(C[0, 0])]
            for t, C in zip(times_corr, integrate(autocorrelation_rhs, covariances[-1], times_corr))]
    return xj, z1z1


# ----------------------------------------------------------
# general case: full density matrix

def dense_solution(N, J, Jxx, Jzzz, hx, gamma, times, times_corr):
    """Same output as free_fermion_solution, from the Lindblad equation of the full density matrix."""
    X, Z = np.array([[0., 1.], [1., 0.]]), np.diag([1., -1.])

    def site(o, j):
        return reduce(np.kron, [o if k == j else np.eye(2) for k in range(N)])

    H = sum(-J * site(Z, j) @ site(Z, j + 1) + Jxx * site(X, j) @ site(X, j + 1) for j in range(N - 1))
    H = H + sum(Jzzz * site(Z, j) @ site(Z, j + 2) for j in range(N - 2))
    H = H + sum(hx * site(X, j) for j in range(N))

    _, V = np.linalg.eigh(H)
    rho0 = np.outer(V[:, 0], V[:, 0].conj())
    jumps = [np.sqrt(gamma) * site(Z, 0)]

    rhos = lindblad_evolution(rho0, H, jumps, times)
    xj = [[t] + [np.real(np.trace(site(X, j) @ rho) / np.trace(rho)) for j in range(N)] for t, rho in zip(times, rhos)]

    rho_ness = rhos[-1] / np.trace(rhos[-1])
    Z1 = site(Z, 0)
    z1z1 = []
    for t, A in zip(times_corr, lindblad_evolution(Z1 @ rho_ness, H, jumps, times_corr)):
        c = np.trace(Z1 @ A)
        z1z1.append([t, c.real, c.imag, abs(c)])
    return xj, z1z1


def main(input_file):
    p = read_input(input_file)
    N, hx, gamma = p.get('N', 10), p.get('hx', 0.5), p.get('gamma', 0.5)
    Jxx, Jzzz = p.get('Jxx', 0.2), p.get('Jzzz', 0.)
    Tness, T, dt, max_dim = p.get('Tness', 5.), p.get('T', 5.), p.get('dt', 0.05), p.get('max_dim', 128)
    J = 1.
    t_measure, t_corr = p.get('t_measure', 0.2), p.get('t_corr', 0.05)   # as in the TN program

    # measurement times of the TN program, plus Tness (start of the autocorrelation)
    times = np.arange(t_measure, Tness + 1e-9, t_measure)
    times_with_ness = np.r_[times, Tness] if times[-1] < Tness - 1e-9 else times
    times_corr = np.arange(t_corr, T + 1e-9, t_corr)

    if Jxx == 0 and Jzzz == 0:
        xj, z1z1 = free_fermion_solution(N, J, hx, gamma, times_with_ness, times_corr)
    elif N <= MAX_DENSE_SITES:
        xj, z1z1 = dense_solution(N, J, Jxx, Jzzz, hx, gamma, times_with_ness, times_corr)
    else:
        sys.exit(f'Jxx or Jzzz != 0 needs the full density matrix: use N <= {MAX_DENSE_SITES}.')
    xj = xj[:len(times)]

    run_dir = run_directory('impurity_N%d_Jxx%.3f_Jzzz%.3f_hx%.3f_hz%.3f_gamma%.3f_Tness%g_T%g_dt%.4f_D%d' % (
        N, Jxx, Jzzz, hx, p.get('hz', 0.), gamma, Tness, T, dt, max_dim))
    return {
        'xj':   save_and_compare(xj, 't ' + ' '.join(f'X_{j}' for j in range(1, N + 1)), run_dir + 'xj.txt'),
        'z1z1': save_and_compare(z1z1, 't Re Im abs', run_dir + 'z1z1.txt'),
    }


if __name__ == '__main__':
    main(sys.argv[1])
