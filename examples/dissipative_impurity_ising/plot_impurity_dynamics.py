"""Plot the output of impurity_dynamics (see ../plot_utils.py): bond dimension and Tr(rho), profile <X_j> during the relaxation, autocorrelation <Z_1(t) Z_1(0)>.

Usage: python3 plot_impurity_dynamics.py [file prefix of the run] [--show]
"""
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import plot_utils as pu  # noqa: E402

import glob  # noqa: E402

prefix, show = pu.parse_arguments('data/impurity_*_xj.txt', '_xj.txt')
ed = prefix + '_exact_diagonalization'
name = os.path.basename(prefix)
pu.save(pu.plot_columns(prefix + '.txt', title=name), prefix, 'observables', show)
pu.save(pu.plot_site_map(prefix + '_xj.txt', '<X_j>', ed + '_xj.txt', title=name), prefix, 'xj', show)
for path in glob.glob(prefix + '_Tness*_z1z1.txt'):
    pu.save(pu.plot_columns(path, ed + '_z1z1.txt', title=name), path[:-len('.txt')], 'plot', show)
pu.show_all(show)
