"""Exact diagonalization of the Rydberg chain of rydberg_chain_tebd, to check the TN results.

    H = sum_j Omega S^x_j + sum_j Delta_j n_j + sum_j V_j n_j n_{j+1} + sum_j V_{j,j+2} n_j n_{j+2},
    n = |down_z><down_z|, Delta_j = -V1 (j even), -V2 (j odd), V_{j,j+2} = 1/(V_j^(-1/6) + V_{j+1}^(-1/6))^6,
starting from the kink |1...1 0...0> with M excitations.
With position disorder (sigmax > 0) the couplings V_j are read from the TN output <root>_Vj.txt.

Usage: python3 rydberg_chain_exact_diagonalization.py input_rydberg_chain_tebd.txt
Output: data/<root>_exact_diagonalization.txt (t, fidelity, half-chain entropy, n_1, ..., n_N),
        compared with the TN files <root>.txt and <root>_nj.txt.
"""
import sys
import numpy as np
import quimb as qu
sys.path.insert(0, '..')
from exact_diagonalization_tools import read_input, save_and_compare

p = read_input(sys.argv[1])
N, M = p.get('N', 12), p.get('M', 2)
V2, Omega = p.get('V2', 2.), p.get('Omega', 0.1)
T, dt, maxDim = p.get('T', 10.), p.get('dt', 0.05), p.get('maxDim', 64)
sigmax, seed = p.get('sigmax', 0.), p.get('seed', 1)
V1 = 1.
t_measure = 10 * dt     # the TN program measures every 10 steps

root = 'data/rydberg_N%d_M%d_V2_%.2f_Om_%.3f_D%d' % (N, M, V2, Omega, maxDim)
if sigmax > 0:
    root += '_sigmax%.5f_seed%d' % (sigmax, seed)
    Vj = np.loadtxt(root + '_Vj.txt')    # disorder realization of the TN run
else:
    d1, d2 = V1 ** (-1 / 6), V2 ** (-1 / 6)
    Vj = np.array([1 / (d1 if j % 2 == 0 else d2) ** 6 for j in range(N - 1)])

# Hamiltonian (site j = 1..N is index j-1)
dims = [2] * N
n = qu.qu(np.diag([0., 1.]), qtype='dop')
Sx = qu.pauli('X') / 2
H = 0
for j in range(1, N + 1):
    Delta = -V1 if j % 2 == 0 else -V2
    H = H + Omega * qu.ikron(Sx, dims, j - 1) + Delta * qu.ikron(n, dims, j - 1)
for j in range(N - 1):
    H = H + Vj[j] * qu.ikron([n, n], dims, [j, j + 1])
for j in range(N - 2):
    VNNN = 1 / (Vj[j] ** (-1 / 6) + Vj[j + 1] ** (-1 / 6)) ** 6
    H = H + VNNN * qu.ikron([n, n], dims, [j, j + 2])

# kink |1...1 0...0>, |1> = |down_z>
psi0 = qu.kron(*[qu.basis_vec(1 if j < M else 0, 2) for j in range(N)])

ts = np.arange(0, T + 1e-9, t_measure)
evo = qu.Evolution(psi0, H, method='solve')
data = []
for t, psi in zip(ts, evo.at_times(ts)):
    fidelity = abs(qu.fidelity(psi, psi0)) ** 2   # quimb returns |<psi|psi0>|
    S = qu.entropy_subsys(psi, dims, sysa=range(N // 2))
    nj = [np.real(qu.expec(qu.ikron(n, dims, j), psi)) for j in range(N)]
    data.append([t, fidelity, S] + nj)
data = np.array(data)

header = 't fidelity S_half ' + ' '.join(f'n_{j}' for j in range(1, N + 1))
save_and_compare(root + '_exact_diagonalization_obs.txt', data[:, :3], 't fidelity S_half', root + '.txt', [1, 2])
save_and_compare(root + '_exact_diagonalization_nj.txt', data[:, [0] + list(range(3, 3 + N))],
                 't ' + ' '.join(f'n_{j}' for j in range(1, N + 1)), root + '_nj.txt')
