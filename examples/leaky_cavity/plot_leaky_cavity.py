"""Plot the output of leaky_cavity (see ../plot_utils.py): Tr(rho), first-spin magnetizations, photon number, bond dimension; profiles <X_j> and <Z_j>.

Usage: python3 plot_leaky_cavity.py [file prefix of the run] [--show]
"""
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import plot_utils as pu  # noqa: E402

prefix, show = pu.parse_arguments('data/leaky_cavity_*_xj.txt', '_xj.txt')
name = os.path.basename(prefix)
pu.save(pu.plot_columns(prefix + '_obs.txt', prefix + '_exact_diagonalization_obs.txt', [2, 3, 4], title=name), prefix, 'observables', show)
pu.save(pu.plot_site_map(prefix + '_xj.txt', '<X_j>', title=name), prefix, 'xj', show)
pu.save(pu.plot_site_map(prefix + '_zj.txt', '<Z_j>', title=name), prefix, 'zj', show)
pu.show_all(show)
