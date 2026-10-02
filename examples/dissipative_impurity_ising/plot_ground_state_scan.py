"""Plot the output of ground_state_scan (see ../plot_utils.py): energy, variance and central magnetization against hx.

Usage: python3 plot_ground_state_scan.py [file prefix of the run] [--show]
"""
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import plot_utils as pu  # noqa: E402

prefix, show = pu.parse_arguments('data/ground_state_scan_*.txt', '.txt')
pu.save(pu.plot_columns(prefix + '.txt', title=os.path.basename(prefix), marker='o'), prefix, 'scan', show)
pu.show_all(show)
