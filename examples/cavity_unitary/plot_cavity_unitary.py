"""Plot the output of cavity_unitary (see ../plot_utils.py): fidelity, collective spin, photon number and bond dimension.

Usage: python3 plot_cavity_unitary.py [run name] [--show]      (default: the latest run in data/)
"""
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import plot_utils as pu  # noqa: E402

data, plots, show = pu.parse_arguments('cavity_unitary_*')
pu.save(pu.plot_columns(data + 'observables.txt', title=os.path.basename(data[:-1])), plots, 'observables', show)
pu.show_all(show)
