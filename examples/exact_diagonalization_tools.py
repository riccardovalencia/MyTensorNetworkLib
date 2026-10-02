"""Shared helpers for the exact-diagonalization (ED) scripts of the examples.

The ED scripts (<example>/<name>_exact_diagonalization.py) read the same input files as the
tensor-network (TN) programs, solve the same problem exactly for small systems and compare the
results with the TN output, if present. Each run writes into data/<run>/ (run name encoding the
parameters); the ED result of a TN file X.txt is written next to it as X_exact_diagonalization.txt. Each script defines main(input_file), which returns
the comparisons ({output: {column: max |TN - ED|}}), so that it can be used by the tests (tests/).
"""
import os
import re
import numpy as np


def read_input(path):
    """Read an ITensor InputGroup file ("input { key = value ... }") into a dict.

    Values are converted to int or float when possible, otherwise kept as strings; integers with
    leading zeros (e.g. a configuration "0110") stay strings.
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
        if re.fullmatch(r'0\d+', value):
            params[key] = value
            continue
        for convert in (int, float):
            try:
                value = convert(value)
                break
            except ValueError:
                pass
        params[key] = value
    return params


def run_directory(run_name):
    """Folder data/<run_name>/ of a run (as built by make_run_directory in the TN programs)."""
    return os.path.join('data', run_name) + '/'


def ed_file_name(tn_file):
    """ED counterpart X_exact_diagonalization.txt of a TN file X.txt (same folder)."""
    base, extension = os.path.splitext(tn_file)
    return f'{base}_exact_diagonalization{extension}'


def compare_tables(ed, tn, columns=None):
    """Maximum absolute difference between ED and TN tables, column by column.

    ed, tn: arrays with one row per time (or theta), first column = time (or theta).
    columns: indices of the TN columns to compare with the ED columns 1, 2, ... (default: same order).
    Rows are matched by their first column; returns (number of common rows, list of differences),
    with an empty list if there are no common rows.
    """
    ed, tn = np.asarray(ed, dtype=float), np.asarray(tn, dtype=float)
    columns = columns if columns is not None else list(range(1, ed.shape[1]))
    rows = [(i, np.argmin(np.abs(tn[:, 0] - t))) for i, t in enumerate(ed[:, 0])]
    rows = [(i, k) for i, k in rows if abs(tn[k, 0] - ed[i, 0]) < 1e-8]
    if not rows:
        return 0, []
    i_ed, i_tn = map(list, zip(*rows))
    return len(rows), [np.abs(tn[i_tn, c_tn] - ed[i_ed, c_ed]).max() for c_ed, c_tn in zip(range(1, ed.shape[1]), columns)]


def save_and_compare(data, header, tn_file, columns=None):
    """Save the ED data next to the TN file (ed_file_name) and compare it with tn_file, if it exists.

    data, TN data: one row per time, first column = time.
    columns: indices of the TN columns to compare with the ED columns 1, 2, ... (default: same order).
    Prints and returns {column name: max |TN - ED|} over the rows with the same first column;
    returns None if tn_file does not exist or has no row in common with the ED data.
    """
    ed_file = ed_file_name(tn_file)
    os.makedirs(os.path.dirname(ed_file) or '.', exist_ok=True)
    np.savetxt(ed_file, data, header=header)
    print(f'ED results written to {ed_file}')
    if not os.path.exists(tn_file):
        print(f'TN results not found ({tn_file}): run the TN program first to compare.')
        return None
    number_rows, differences = compare_tables(data, np.loadtxt(tn_file, ndmin=2), columns)
    if number_rows == 0:
        print('No common rows (same first column) between ED and TN results.')
        return None
    names = header.split()
    names = [names[c] if c < len(names) else f'column {c}' for c in range(1, len(differences) + 1)]
    print(f'Comparison with {tn_file} ({number_rows} common rows): max |TN - ED|')
    for name, diff in zip(names, differences):
        print(f'  {name:>12s}: {diff:.2e}')
    return dict(zip(names, differences))


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
