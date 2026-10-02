"""Plot the output of rydberg_chain_tebd (see ../plot_utils.py): fidelity, half-chain entropy and bond dimension; Rydberg densities n_j and entropy of each cut S_j (heatmaps, profiles, time traces).

Usage: python3 plot_rydberg_chain_tebd.py [run name] [--show]      (default: the latest run in data/)
"""
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import plot_utils as pu  # noqa: E402

data, plots, show = pu.parse_arguments('rydberg_*')
run = os.path.basename(data[:-1])
pu.save(pu.plot_columns(data + 'observables.txt', title=run), plots, 'observables', show)
pu.save(pu.plot_site_map(data + 'nj.txt', 'n_j', title=run), plots, 'nj', show)
pu.save(pu.plot_site_map(data + 'entropy.txt', 'S_j', site_label='cut (j, j+1)', title=run), plots, 'entropy', show)
pu.show_all(show)
