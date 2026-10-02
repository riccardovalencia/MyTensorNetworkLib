"""Plot the output of cavity_unitary (see ../plot_utils.py): fidelity, collective spin, photon number and bond dimension.

Usage: python3 plot_cavity_unitary.py [file prefix of the run] [--show]
"""
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import plot_utils as pu  # noqa: E402

prefix, show = pu.parse_arguments('data/cavity_unitary_*.txt', '.txt')
pu.save(pu.plot_columns(prefix + '.txt', prefix + '_exact_diagonalization.txt', title=os.path.basename(prefix)), prefix, 'observables', show)
pu.show_all(show)
