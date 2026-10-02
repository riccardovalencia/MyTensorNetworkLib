"""Plot the output of rydberg_chain_tebd (see ../plot_utils.py): fidelity, half-chain entropy and bond dimension; Rydberg densities n_j and entropy of each cut S_j (heatmaps, profiles, time traces).

Usage: python3 plot_rydberg_chain_tebd.py [file prefix of the run] [--show]
"""
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import plot_utils as pu  # noqa: E402

prefix, show = pu.parse_arguments('data/rydberg_*_nj.txt', '_nj.txt')
ed = prefix + '_exact_diagonalization'
name = os.path.basename(prefix)
pu.save(pu.plot_columns(prefix + '.txt', ed + '_obs.txt', title=name), prefix, 'observables', show)
pu.save(pu.plot_site_map(prefix + '_nj.txt', 'n_j', ed + '_nj.txt', title=name), prefix, 'nj', show)
pu.save(pu.plot_site_map(prefix + '_entropy.txt', 'S_j', ed + '_entropy.txt', site_label='cut (j, j+1)', title=name), prefix, 'entropy', show)
pu.show_all(show)
