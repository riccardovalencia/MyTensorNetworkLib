"""Plot the output of spin_chain_ground_state (see ../plot_utils.py): profiles <Z_j>, <Z_j Z_{j+1}> and entanglement entropy of each cut, energy of every DMRG restart.

Usage: python3 plot_spin_chain_ground_state.py [run name] [--show]      (default: the latest run in data/)
"""
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import plot_utils as pu  # noqa: E402

data, plots, show = pu.parse_arguments('spin_chain_*')
energy, variance, converged = pu.load(data + 'energy.txt')[1][0]
title = f'{os.path.basename(data[:-1])}\nE = {energy:.10f}, variance = {variance:.1e}' + ('' if converged else ' (not converged)')
pu.save(pu.plot_columns(data + 'profile.txt', title=title, marker='o'), plots, 'profile', show)
pu.save(pu.plot_columns(data + 'restarts.txt', title=title, marker='o'), plots, 'restarts', show)
pu.show_all(show)
