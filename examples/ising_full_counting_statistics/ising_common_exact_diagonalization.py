"""Dense exact diagonalization of the Ising chain of ising_quench / ising_thermal (small N).

    H = -J sum_j [ X_j X_{j+1} + hx X_j + hz Z_j ]
"""
from functools import reduce
import numpy as np

X = np.array([[0., 1.], [1., 0.]])
Z = np.diag([1., -1.])
HADAMARD = np.array([[1., 1.], [1., -1.]]) / np.sqrt(2)


def site_operator(o, j, N):
    """Operator o on site j (0-indexed) of an N-site chain."""
    return reduce(np.kron, [o if k == j else np.eye(2) for k in range(N)])


def ising_hamiltonian(N, J, hx, hz):
    H = sum(-J * site_operator(X, j, N) @ site_operator(X, j + 1, N) for j in range(N - 1))
    return H + sum(-J * hx * site_operator(X, j, N) - J * hz * site_operator(Z, j, N) for j in range(N))


def block_start(N, l):
    """First site (1-indexed) of the block of l sites, centered as in the TN code."""
    size = l - 1
    return N // 2 - size // 2 if size % 2 == 0 else N // 2 - (size + 1) // 2


def apply_hadamard(tensor, axes):
    """Apply the Hadamard matrix on the given axes of a tensor with legs of dimension 2."""
    for axis in axes:
        tensor = np.moveaxis(np.tensordot(HADAMARD, tensor, axes=(1, axis)), 0, axis)
    return tensor


def x_basis_probabilities(rho_or_psi, N):
    """Probabilities of the product states of the x basis for a pure state (vector) or a density
    matrix, as an array of shape (2,)*N (index 0 = |+x>, 1 = |-x>)."""
    if rho_or_psi.ndim == 1:
        return np.abs(apply_hadamard(rho_or_psi.reshape((2,) * N), range(N))) ** 2
    rho_x = apply_hadamard(rho_or_psi.reshape((2,) * (2 * N)), range(2 * N))
    d = 2 ** N
    return np.real(np.diag(rho_x.reshape(d, d))).reshape((2,) * N)


def generating_function(probabilities, N, l, theta):
    """G(theta) = <exp(i theta S^x_A)>, S^x_A = sum_{j in A} X_j / 2, from the x-basis probabilities."""
    start = block_start(N, l) - 1
    p = probabilities.sum(axis=tuple(k for k in range(N) if not start <= k < start + l))
    minus_count = np.add.reduce(np.indices(p.shape), axis=0)   # number of |-x> spins in the block
    sx = l / 2 - minus_count
    return np.array([np.sum(p * np.exp(1j * t * sx)) for t in theta])


def generating_function_table(rho_or_psi, N, theta):
    """Rows theta, Re G_1, Im G_1, ..., Re G_{N/2}, Im G_{N/2} (format of the TN files)."""
    probabilities = x_basis_probabilities(rho_or_psi, N)
    columns = [theta]
    for l in range(1, N // 2 + 1):
        G = generating_function(probabilities, N, l, theta)
        columns += [G.real, G.imag]
    return np.column_stack(columns)
