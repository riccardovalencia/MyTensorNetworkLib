"""Shared helpers for the exact-diagonalization (ED) scripts of the examples.

The ED scripts (<example>/<name>_exact_diagonalization.py) read the same input files as the
tensor-network (TN) programs, solve the same problem exactly for small systems and compare the
results with the TN output, if present in data/.
"""
import os
import re
import numpy as np


def read_input(path):
    """Read an ITensor InputGroup file ("input { key = value ... }") into a dict.

    Values are converted to int or float when possible, otherwise kept as strings.
    """
    with open(path) as f:
        text = f.read()
    body = text[text.index('{') + 1:text.rindex('}')]
    params = {}
    for line in body.splitlines():
        line = line.split('#')[0].split('//')[0].strip()
        if not line:
            continue
        key, value = [s.strip() for s in line.split('=', 1)]
        for convert in (int, float):
            try:
                value = convert(value)
                break
            except ValueError:
                pass
        params[key] = value
    return params


def tn_root(fmt, *args):
    """C-style file prefix, as built by the TN programs with tinyformat::format."""
    return fmt % args


def save_and_compare(ed_file, data, header, tn_file, columns=None):
    """Save the ED data and compare it with the TN output tn_file, if it exists.

    data, TN data: one row per time, first column = time.
    columns: indices of the TN columns to compare with the ED columns 1, 2, ... (default: same order).
    Rows are matched by their first column (time or theta); the maximum absolute difference of each
    column is printed.
    """
    os.makedirs(os.path.dirname(ed_file) or '.', exist_ok=True)
    np.savetxt(ed_file, data, header=header)
    print(f'ED results written to {ed_file}')
    if not os.path.exists(tn_file):
        print(f'TN results not found ({tn_file}): run the TN program first to compare.')
        return
    tn = np.loadtxt(tn_file, ndmin=2)
    ed = np.asarray(data)
    columns = columns if columns is not None else list(range(1, ed.shape[1]))
    names = header.split()
    rows = [(i, np.argmin(np.abs(tn[:, 0] - t))) for i, t in enumerate(ed[:, 0])]
    rows = [(i, k) for i, k in rows if abs(tn[k, 0] - ed[i, 0]) < 1e-8]
    if not rows:
        print('No common rows (same first column) between ED and TN results.')
        return
    i_ed, i_tn = zip(*rows)
    print(f'Comparison with {tn_file} ({len(rows)} common rows): max |TN - ED|')
    for c_ed, c_tn in zip(range(1, ed.shape[1]), columns):
        diff = np.abs(tn[list(i_tn), c_tn] - ed[list(i_ed), c_ed]).max()
        name = names[c_ed] if c_ed < len(names) else f'column {c_ed}'
        print(f'  {name:>12s}: {diff:.2e}')


def lindblad_evolution(rho0, H, jump_operators, times, rtol=1e-10, atol=1e-12):
    """Density matrices rho(t) at the given times for d rho/dt = -i[H, rho] + sum_k D[L_k] rho,
    D[L] rho = L rho L^dag - 1/2 {L^dag L, rho} (include the rates in the jump operators).

    Integrates the vectorized (row-major) Lindblad equation with scipy's solve_ivp.
    """
    import scipy.sparse as sp
    from scipy.integrate import solve_ivp
    d = H.shape[0]
    I = sp.identity(d, format='csr')
    H = sp.csr_matrix(H)
    # vec(A rho B) = (A kron B^T) vec(rho) for row-major vectorization
    L = -1j * (sp.kron(H, I) - sp.kron(I, H.T))
    for Lk in jump_operators:
        Lk = sp.csr_matrix(Lk)
        LdL = (Lk.conj().T @ Lk).tocsr()
        L = L + sp.kron(Lk, Lk.conj()) - 0.5 * (sp.kron(LdL, I) + sp.kron(I, LdL.T))
    L = L.tocsr()
    sol = solve_ivp(lambda t, y: L @ y, (0, times[-1]), np.asarray(rho0, dtype=complex).ravel(),
                    t_eval=times, method='DOP853', rtol=rtol, atol=atol)
    return [sol.y[:, k].reshape(d, d) for k in range(len(times))]


def theta_grid(number_points):
    """Values of theta used by the TN generating functions (analysis/full_counting_statistics.h)."""
    theta = [-np.pi]
    for col in range(number_points - 1):
        if col <= number_points // 4 or col >= 3 * number_points // 4 - 1:
            theta.append(theta[-1] + (np.pi - 1.) / (number_points // 4))
        else:
            theta.append(theta[-1] + 2. / (number_points / 2.))
    return np.array(theta)
