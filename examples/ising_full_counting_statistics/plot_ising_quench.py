"""Plot the output of ising_quench (see ../plot_utils.py): magnetization <X_j>(t) of each spin and
entanglement entropy of each cut (heatmaps, profiles with the thermal <X_j>, time traces), distance of the generating functions from the thermal ones against time, and the
generating functions at each time with the thermal ones (red dotted).

Usage: python3 plot_ising_quench.py [run name] [--show]      (default: the latest run in data/)
"""
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import plot_utils as pu  # noqa: E402

data, plots, show = pu.parse_arguments('ising_quench_*')
run = os.path.basename(data[:-1])
pu.save(pu.plot_site_map(data + 'xj.txt', '<X_j>', title=run, reference=data + 'thermal_xj.txt', reference_label='thermal'), plots, 'xj', show)
pu.save(pu.plot_site_map(data + 'entropy.txt', 'S_j', site_label='cut (j, j+1)', title=run), plots, 'entropy', show)
pu.save(pu.plot_columns(data + 'distance_to_thermal.txt', title=run + '\nmax_theta |G_l(theta, t) - G_l^thermal(theta)|', marker='o'),
        plots, 'distance_to_thermal', show)
pu.save(pu.plot_generating_function(data + 'thermal_gf.txt', title=run + ', thermal state'), plots, 'thermal_gf', show)
for path in pu.tn_files(data, 'gf_t*.txt'):
    name = os.path.basename(path)[:-len('.txt')]
    pu.save(pu.plot_generating_function(path, title=f'{run}, {name}', reference=data + 'thermal_gf.txt', reference_label='thermal'),
            plots, name, show)
pu.show_all(show)
