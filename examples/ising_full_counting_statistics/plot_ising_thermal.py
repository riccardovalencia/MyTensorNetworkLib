"""Plot the output of ising_thermal (see ../plot_utils.py): energy density along the imaginary-time evolution and generating functions of the thermal state.

Usage: python3 plot_ising_thermal.py [run name] [--show]      (default: the latest run in data/)
"""
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import plot_utils as pu  # noqa: E402

data, plots, show = pu.parse_arguments('ising_thermal_*')
run = os.path.basename(data[:-1])
pu.save(pu.plot_columns(data + 'energy.txt', title=run), plots, 'energy', show)
pu.save(pu.plot_generating_function(data + 'gf.txt', title=run), plots, 'gf', show)
pu.show_all(show)
