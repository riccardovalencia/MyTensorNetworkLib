"""Plot the output of ising_quench (see ../plot_utils.py): entanglement entropy of each cut (heatmap, profiles, time traces) and generating functions at each time.

Usage: python3 plot_ising_quench.py [run name] [--show]      (default: the latest run in data/)
"""
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import plot_utils as pu  # noqa: E402

data, plots, show = pu.parse_arguments('ising_quench_*')
run = os.path.basename(data[:-1])
pu.save(pu.plot_site_map(data + 'entropy.txt', 'S_j', site_label='cut (j, j+1)', title=run), plots, 'entropy', show)
for path in pu.tn_files(data, 'gf_t*.txt'):
    name = os.path.basename(path)[:-len('.txt')]
    pu.save(pu.plot_generating_function(path, title=f'{run}, {name}'), plots, name, show)
pu.show_all(show)
