"""Plot the output of ground_state_scan (see ../plot_utils.py): energy, variance and central magnetization against hx.

Usage: python3 plot_ground_state_scan.py [run name] [--show]      (default: the latest run in data/)
"""
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import plot_utils as pu  # noqa: E402

data, plots, show = pu.parse_arguments('ground_state_scan_*')
pu.save(pu.plot_columns(data + 'scan.txt', title=os.path.basename(data[:-1]), marker='o'), plots, 'scan', show)
pu.show_all(show)
