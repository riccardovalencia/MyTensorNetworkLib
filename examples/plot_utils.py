"""Shared helpers of the scripts <example>/plot_<program>.py, which plot the output of a TN program.

Each run of a program writes its data into data/<run>/ (run name encoding the parameters); the plot
scripts save the figures of the run into plots/<run>/<name>.png. Without arguments a script plots the
latest run of its program, otherwise the run given by name (or by its data folder); --show also opens
the figures. Exact-diagonalization (ED) results are optional: the ED counterpart X_exact_diagonalization.txt
of a TN file X.txt is overlaid as dashed lines when it exists and matches the TN data, and skipped otherwise.

    python3 plot_rydberg_chain_tebd.py                                    # latest run
    python3 plot_rydberg_chain_tebd.py rydberg_N12_M2_V2_2.00_Om_0.300_D64 --show
"""
import argparse
import glob
import os
import numpy as np
import matplotlib
from exact_diagonalization_tools import ed_file_name


def parse_arguments(pattern):
    """Data folder data/<run>/ and plot folder plots/<run>/ (created) of the run to plot, and --show.

    pattern: glob of the run names of the program (e.g. 'rydberg_*'); without a run argument the run
    with the most recently written file is used.
    """
    parser = argparse.ArgumentParser(description='Plot the output of a TN run.')
    parser.add_argument('run', nargs='?', help='run name, or its folder data/<run> (default: the latest run)')
    parser.add_argument('--show', action='store_true', help='open the figures instead of only saving them')
    args = parser.parse_args()
    if not args.show:
        matplotlib.use('Agg')
    if args.run:
        run = os.path.basename(os.path.normpath(args.run))
    else:
        runs = [d for d in glob.glob(os.path.join('data', pattern)) if os.path.isdir(d) and os.listdir(d)]
        if not runs:
            raise SystemExit(f'No run data/{pattern}: run the program first.')
        run = os.path.basename(max(runs, key=lambda d: max(os.path.getmtime(f) for f in glob.glob(d + '/*'))))
    data_dir = os.path.join('data', run) + '/'
    if not os.path.isdir(data_dir):
        raise SystemExit(f'No folder {data_dir}.')
    plots_dir = os.path.join('plots', run) + '/'
    os.makedirs(plots_dir, exist_ok=True)
    return data_dir, plots_dir, args.show


def load(path):
    """Column names (from the header "# a . b . c") and the data (2D array) of a TN/ED output file."""
    with open(path) as f:
        header = f.readline()
    names = [n.strip() for n in header.lstrip('#').split(' . ')] if header.startswith('#') else []
    return names, np.loadtxt(path, ndmin=2)


def load_exact_diagonalization(path, number_columns=None):
    """Data of the ED counterpart of the TN file path, or None if it does not exist or has a different
    number of columns than the TN data (e.g. a file left by an ED run with other parameters)."""
    ed = ed_file_name(path)
    if not os.path.exists(ed):
        return None
    data = load(ed)[1]
    if number_columns is not None and data.shape[1] != number_columns:
        print(f'skipped {ed}: {data.shape[1]} columns instead of {number_columns}')
        return None
    return data


def plot_columns(path, ed_columns=None, title='', marker=''):
    """One panel per column against the first one (time, field, site, ...), with the ED result if any.

    ed_columns: TN columns matching the ED columns 1, 2, ... (default: the same order).
    """
    import matplotlib.pyplot as plt
    names, data = load(path)
    columns = range(1, data.shape[1])
    ncols = min(3, len(columns))
    nrows = -(-len(columns) // ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.5 * ncols, 3.2 * nrows), squeeze=False)
    ed_data = load_exact_diagonalization(path)
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


def plot_site_map(path, quantity, site_label='site j', title=''):
    """Site-resolved quantity in a file with rows "t q_1 q_2 ...": heatmap over (site, t), profiles
    at a few times and the time trace of every site (ED dashed, if any)."""
    import matplotlib.pyplot as plt
    data = load(path)[1]
    t, values = data[:, 0], data[:, 1:]
    sites = np.arange(1, values.shape[1] + 1)
    ed_data = load_exact_diagonalization(path, data.shape[1])

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


def plot_generating_function(path, title=''):
    """Generating functions G_l(theta) of a file with rows "theta Re G_1 Im G_1 Re G_2 ...":
    real and imaginary parts, one line per block size l (ED dashed, if any)."""
    import matplotlib.pyplot as plt
    data = load(path)[1]
    ed_data = load_exact_diagonalization(path, data.shape[1])
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


def tn_files(data_dir, pattern):
    """TN files of data_dir matching pattern (e.g. 'gf_t*.txt'), without their ED counterparts, sorted."""
    return sorted(f for f in glob.glob(data_dir + pattern) if not f.endswith('_exact_diagonalization.txt'))


def save(fig, plots_dir, name, show):
    """Save the figure as plots_dir/<name>.png; with show, keep it open for show_all."""
    path = f'{plots_dir}{name}.png'
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
