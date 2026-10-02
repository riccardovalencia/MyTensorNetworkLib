"""Shared helpers of the scripts <example>/plot_<program>.py, which plot the output of a TN program.

Each script plots the most recent run found in data/ (or the run whose file prefix, the part of the
data file names before the suffix, is given as argument), saves the figures as <prefix>_<name>.png
and, with --show, opens them. Exact-diagonalization results, if present, are overlaid as dashed lines.

    python3 plot_rydberg_chain_tebd.py                              # latest run
    python3 plot_rydberg_chain_tebd.py data/rydberg_N12_..._D64 --show
"""
import argparse
import glob
import os
import numpy as np
import matplotlib


def parse_arguments(pattern, suffix):
    """Prefix of the run to plot and the --show flag.

    pattern: glob of a data file written by every run (e.g. 'data/rydberg_*_nj.txt'), suffix: its ending
    after the prefix (e.g. '_nj.txt'). Without a prefix argument the newest matching file is used
    (files written by the ED scripts are skipped).
    """
    parser = argparse.ArgumentParser(description='Plot the output of a TN run.')
    parser.add_argument('prefix', nargs='?', help='file prefix of the run (default: the latest run in data/)')
    parser.add_argument('--show', action='store_true', help='open the figures instead of only saving them')
    args = parser.parse_args()
    if not args.show:
        matplotlib.use('Agg')
    if args.prefix:
        return args.prefix, args.show
    files = [f for f in glob.glob(pattern) if 'exact_diagonalization' not in f and not f.endswith('_ED.txt')]
    if not files:
        raise SystemExit(f'No file {pattern}: run the program first.')
    return max(files, key=os.path.getmtime)[:-len(suffix)], args.show


def load(path):
    """Column names (from the header "# a . b . c") and the data (2D array) of a TN/ED output file."""
    with open(path) as f:
        header = f.readline()
    names = [n.strip() for n in header.lstrip('#').split(' . ')] if header.startswith('#') else []
    return names, np.loadtxt(path, ndmin=2)


def existing(path):
    """path if the file exists, otherwise None (for optional ED files)."""
    return path if path and os.path.exists(path) else None


def plot_columns(path, ed=None, ed_columns=None, title='', marker=''):
    """One panel per column against the first one (time, field, site, ...).

    ed: optional ED file with the same first column; ed_columns: TN columns matching the ED columns
    1, 2, ... (default: the same order).
    """
    import matplotlib.pyplot as plt
    names, data = load(path)
    columns = range(1, data.shape[1])
    ncols = min(3, len(columns))
    nrows = -(-len(columns) // ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.5 * ncols, 3.2 * nrows), squeeze=False)
    ed_data = load(ed)[1] if existing(ed) else None
    ed_of = dict(zip(ed_columns or columns, range(1, ed_data.shape[1]))) if ed_data is not None else {}
    for ax, c in zip(axes.flat, columns):
        ax.plot(data[:, 0], data[:, c], marker=marker, label='TN')
        if c in ed_of:
            ax.plot(ed_data[:, 0], ed_data[:, ed_of[c]], 'k--', lw=1, label='ED')
            ax.legend()
        ax.set_xlabel(names[0] if names else '')
        ax.set_ylabel(names[c] if c < len(names) else f'column {c}')
    for ax in axes.flat[len(columns):]:
        ax.set_visible(False)
    fig.suptitle(title)
    fig.tight_layout()
    return fig


def plot_site_map(path, quantity, ed=None, site_label='site j', title=''):
    """Site-resolved quantity in a file with rows "t q_1 q_2 ...": heatmap over (site, t), profiles
    at a few times and the time trace of every site (ED dashed, if given)."""
    import matplotlib.pyplot as plt
    data = load(path)[1]
    t, values = data[:, 0], data[:, 1:]
    sites = np.arange(1, values.shape[1] + 1)
    ed_data = load(ed)[1] if existing(ed) else None

    fig, (ax_map, ax_profiles, ax_traces) = plt.subplots(1, 3, figsize=(15, 4))
    mesh = ax_map.pcolormesh(sites, t, values, shading='nearest', cmap='viridis')
    fig.colorbar(mesh, ax=ax_map, label=quantity)
    ax_map.set(xlabel=site_label, ylabel='t', title='heatmap')

    snapshots = np.unique(np.linspace(0, len(t) - 1, min(5, len(t))).astype(int))
    colors = plt.cm.plasma(np.linspace(0, 0.9, len(snapshots)))
    for k, color in zip(snapshots, colors):
        ax_profiles.plot(sites, values[k], 'o-', color=color, label=f't = {t[k]:g}')
        if ed_data is not None:
            row = np.argmin(np.abs(ed_data[:, 0] - t[k]))
            if abs(ed_data[row, 0] - t[k]) < 1e-8:
                ax_profiles.plot(sites, ed_data[row, 1:], 'k--', lw=1)
    ax_profiles.set(xlabel=site_label, ylabel=quantity, title='profiles' + (' (ED dashed)' if ed_data is not None else ''))
    ax_profiles.legend(fontsize='small')

    cmap = plt.cm.viridis
    for j in range(len(sites)):
        ax_traces.plot(t, values[:, j], color=cmap(j / max(1, len(sites) - 1)))
        if ed_data is not None:
            ax_traces.plot(ed_data[:, 0], ed_data[:, j + 1], 'k--', lw=0.6)
    fig.colorbar(plt.cm.ScalarMappable(norm=matplotlib.colors.Normalize(sites[0], sites[-1]), cmap=cmap),
                 ax=ax_traces, label=site_label)
    ax_traces.set(xlabel='t', ylabel=quantity, title='time traces')
    fig.suptitle(title)
    fig.tight_layout()
    return fig


def plot_generating_function(path, ed=None, title=''):
    """Generating functions G_l(theta) of a file with rows "theta Re G_1 Im G_1 Re G_2 ...":
    real and imaginary parts, one line per block size l (ED dashed)."""
    import matplotlib.pyplot as plt
    data = load(path)[1]
    ed_data = load(ed)[1] if existing(ed) else None
    blocks = (data.shape[1] - 1) // 2
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.6))
    colors = plt.cm.viridis(np.linspace(0, 0.9, blocks))
    for part, ax in enumerate(axes):
        for l, color in zip(range(1, blocks + 1), colors):
            ax.plot(data[:, 0], data[:, 2 * l - 1 + part], color=color, label=f'l = {l}')
            if ed_data is not None:
                ax.plot(ed_data[:, 0], ed_data[:, 2 * l - 1 + part], 'k--', lw=0.6)
        ax.set(xlabel='theta', ylabel=('Re' if part == 0 else 'Im') + ' G_l(theta)')
    axes[0].legend(fontsize='small')
    fig.suptitle(title)
    fig.tight_layout()
    return fig


def save(fig, prefix, name, show):
    """Save the figure as <prefix>_<name>.png; with show, keep it open for plt.show()."""
    path = f'{prefix}_{name}.png'
    fig.savefig(path, dpi=120)
    print(f'saved {path}')
    if not show:
        import matplotlib.pyplot as plt
        plt.close(fig)


def show_all(show):
    """Open the figures if --show was given."""
    if show:
        import matplotlib.pyplot as plt
        plt.show()
