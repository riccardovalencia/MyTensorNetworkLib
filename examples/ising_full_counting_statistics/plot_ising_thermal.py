"""Plot the output of ising_thermal (see ../plot_utils.py): energy density along the imaginary-time evolution and generating functions of the thermal state.

Usage: python3 plot_ising_thermal.py [file prefix of the run] [--show]
"""
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import plot_utils as pu  # noqa: E402

prefix, show = pu.parse_arguments('data/ising_thermal_*_energy.txt', '_energy.txt')
name = os.path.basename(prefix)
pu.save(pu.plot_columns(prefix + '_energy.txt', title=name), prefix, 'energy', show)
pu.save(pu.plot_generating_function(prefix + '_gf.txt', prefix + '_exact_diagonalization_gf.txt', title=name), prefix, 'gf', show)
pu.show_all(show)
