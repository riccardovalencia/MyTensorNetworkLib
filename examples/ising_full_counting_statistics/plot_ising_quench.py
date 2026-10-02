"""Plot the output of ising_quench (see ../plot_utils.py): entanglement entropy of each cut (heatmap, profiles, time traces) and generating functions at each time.

Usage: python3 plot_ising_quench.py [file prefix of the run] [--show]
"""
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import plot_utils as pu  # noqa: E402

import glob  # noqa: E402

prefix, show = pu.parse_arguments('data/ising_quench_*_entropy.txt', '_entropy.txt')
ed = prefix + '_exact_diagonalization'
name = os.path.basename(prefix)
pu.save(pu.plot_site_map(prefix + '_entropy.txt', 'S_j', ed + '_entropy.txt', site_label='cut (j, j+1)', title=name), prefix, 'entropy', show)
for path in sorted(glob.glob(prefix + '_gf_t*.txt')):
    time = path[len(prefix) + len('_gf_'):-len('.txt')]
    pu.save(pu.plot_generating_function(path, f'{ed}_gf_{time}.txt', title=f'{name}, {time}'), prefix, f'gf_{time}', show)
pu.show_all(show)
