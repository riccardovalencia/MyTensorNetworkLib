"""Plot the output of leaky_cavity (see ../plot_utils.py): Tr(rho), first-spin magnetizations, photon number, bond dimension; profiles <X_j> and <Z_j>.

Usage: python3 plot_leaky_cavity.py [run name] [--show]      (default: the latest run in data/)
"""
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import plot_utils as pu  # noqa: E402

data, plots, show = pu.parse_arguments('leaky_cavity_*')
run = os.path.basename(data[:-1])
pu.save(pu.plot_columns(data + 'observables.txt', ed_columns=[2, 3, 4], title=run), plots, 'observables', show)
pu.save(pu.plot_site_map(data + 'xj.txt', '<X_j>', title=run), plots, 'xj', show)
pu.save(pu.plot_site_map(data + 'zj.txt', '<Z_j>', title=run), plots, 'zj', show)
pu.show_all(show)
