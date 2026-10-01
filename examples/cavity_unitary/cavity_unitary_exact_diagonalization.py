"""Exact diagonalization of the closed cavity model of cavity_unitary, to check the TN results.

    H = omega0 a^dag a + h sum_j Z_j + (g/sqrt(N)) * coupling,
    coupling = (a + a^dag) sum_j X_j ("dicke") or sum_j (sigma^+_j a + sigma^-_j a^dag) ("tavis"),
with the boson truncated at max_occ and the initial state |0> (x) (cos(theta/2)|up_z> + sin(theta/2)|down_z>)^N.

Usage: python3 cavity_unitary_exact_diagonalization.py input_cavity_unitary.txt
Output: data/<root>_exact_diagonalization.txt (t, fidelity, <S^x>/N, <S^z>/N, <a^dag a>/N),
        compared with the TN file <root>.txt.
"""
import os
import sys
import numpy as np
import quimb as qu
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from exact_diagonalization_tools import read_input, save_and_compare

p = read_input(sys.argv[1])
N, max_occ = p.get('N', 3), p.get('max_occ', 6)
omega0, h, g = p.get('omega0', 1.), p.get('h', 0.5), p.get('g', 1.5)
theta = p.get('theta', 0.5) * np.pi
T, dt, coupling = p.get('T', 50.), p.get('dt', 0.005), p.get('coupling', 'dicke')
t_measure = 10 * dt     # the TN program measures every 10 steps

# operators: site 0 = boson, sites 1..N = spins (|up_z> = basis state 0)
dims = [max_occ + 1] + [2] * N
a = qu.qu(np.diag(np.sqrt(np.arange(1, max_occ + 1)), k=1), qtype='dop')
X, Z = qu.pauli('X'), qu.pauli('Z')
sigma_minus = qu.qu([[0, 0], [1, 0]], qtype='dop')   # |up_z> -> |down_z>

def op(o, j):
    return qu.ikron(o, dims, j)

H = omega0 * op(a.H @ a, 0)
for j in range(1, N + 1):
    H = H + h * op(Z, j)
    if coupling == 'dicke':
        H = H + g / np.sqrt(N) * qu.ikron([a + a.H, X], dims, [0, j])
    else:
        H = H + g / np.sqrt(N) * (qu.ikron([a, sigma_minus.H], dims, [0, j]) + qu.ikron([a.H, sigma_minus], dims, [0, j]))

spin = np.cos(theta / 2) * qu.basis_vec(0, 2) + np.sin(theta / 2) * qu.basis_vec(1, 2)
psi0 = qu.kron(qu.basis_vec(0, max_occ + 1), *[spin] * N)

Sx = sum(op(X, j) for j in range(1, N + 1)) / N
Sz = sum(op(Z, j) for j in range(1, N + 1)) / N
n_photon = op(a.H @ a, 0) / N

times = np.arange(t_measure, T + 1e-9, t_measure)
evo = qu.Evolution(psi0, H, method='solve')
data = [[t, abs(qu.fidelity(psi, psi0)) ** 2] + [np.real(qu.expec(O, psi)) for O in (Sx, Sz, n_photon)]
        for t, psi in zip(times, evo.at_times(times))]

root = 'data/cavity_unitary_%s_N%d_maxocc%d_omega%.2f_h%.2f_g%.2f' % (coupling, N, max_occ, omega0, h, g)
save_and_compare(root + '_exact_diagonalization.txt', data, 't fidelity Sx/N Sz/N n_photon/N', root + '.txt')
