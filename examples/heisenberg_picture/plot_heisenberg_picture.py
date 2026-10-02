"""Plot the output of heisenberg_picture (see ../plot_utils.py): <Z_j(t)> in the Schroedinger and Heisenberg pictures, bond dimensions of the state and of the operator.

Usage: python3 plot_heisenberg_picture.py [file prefix of the run] [--show]
"""
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import plot_utils as pu  # noqa: E402

prefix, show = pu.parse_arguments('data/heisenberg_*.txt', '.txt')
pu.save(pu.plot_columns(prefix + '.txt', prefix + '_ED.txt', [1, 2], title=os.path.basename(prefix)), prefix, 'observables', show)
pu.show_all(show)
