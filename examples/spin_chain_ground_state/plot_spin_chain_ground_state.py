"""Plot the output of spin_chain_ground_state (see ../plot_utils.py): profiles <Z_j>, <Z_j Z_{j+1}> and entanglement entropy of each cut, energy of every DMRG restart.

Usage: python3 plot_spin_chain_ground_state.py [file prefix of the run] [--show]
"""
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import plot_utils as pu  # noqa: E402

prefix, show = pu.parse_arguments('data/spin_chain_*_profile.txt', '_profile.txt')
energy, variance, converged = pu.load(prefix + '.txt')[1][0]
name = f'{os.path.basename(prefix)}\nE = {energy:.10f}, variance = {variance:.1e}' + ('' if converged else ' (not converged)')
pu.save(pu.plot_columns(prefix + '_profile.txt', prefix + '_profile_ED.txt', title=name, marker='o'), prefix, 'profile', show)
pu.save(pu.plot_columns(prefix + '_restarts.txt', title=name, marker='o'), prefix, 'restarts', show)
pu.show_all(show)
