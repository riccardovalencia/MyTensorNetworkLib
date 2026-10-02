"""Plot the output of impurity_dynamics (see ../plot_utils.py): bond dimension and Tr(rho), profile <X_j> during the relaxation, autocorrelation <Z_1(t) Z_1(0)>.

Usage: python3 plot_impurity_dynamics.py [run name] [--show]      (default: the latest run in data/)
"""
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import plot_utils as pu  # noqa: E402

data, plots, show = pu.parse_arguments('impurity_*')
run = os.path.basename(data[:-1])
pu.save(pu.plot_columns(data + 'observables.txt', title=run), plots, 'observables', show)
pu.save(pu.plot_site_map(data + 'xj.txt', '<X_j>', title=run), plots, 'xj', show)
pu.save(pu.plot_columns(data + 'z1z1.txt', title=run), plots, 'z1z1', show)
pu.show_all(show)
