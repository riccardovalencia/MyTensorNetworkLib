"""Plot the output of heisenberg_picture (see ../plot_utils.py): <Z_j(t)> in the Schroedinger and Heisenberg pictures, bond dimensions of the state and of the operator.

Usage: python3 plot_heisenberg_picture.py [run name] [--show]      (default: the latest run in data/)
"""
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import plot_utils as pu  # noqa: E402

data, plots, show = pu.parse_arguments('heisenberg_*')
pu.save(pu.plot_columns(data + 'observables.txt', ed_columns=[1, 2], title=os.path.basename(data[:-1])), plots, 'observables', show)
pu.show_all(show)
