"""Launcher for the pre-registered REPEAT of the pairs run (pairs_repeat/PREREG_PAIRS.md). Heuremen LLC, Wayfinder and Bones, 22 October 2026.
Puts control_run on sys.path (so the copied script's 'import run as base' finds control_run/run.py, unchanged) and executes the BYTE-IDENTICAL
copy pairs_repeat/run_pairs_iq.py with the same argv. The copy's HERE is pairs_repeat, so every output (pairsjob_*, pairsraw_*,
pairs_backend_properties_*) lands here and nothing in control_run is touched; the prereg hash it records is pairs_repeat/PREREG_PAIRS.sha256's first line.
Usage: repeat_pairs.py --chip ibm_fez [--fire]     (dry-run without --fire)"""
import os, sys, runpy
HERE = os.path.dirname(os.path.abspath(__file__)); CTRL = os.path.join(os.path.dirname(HERE), 'control_run'); SCRIPT = os.path.join(HERE, 'run_pairs_iq.py')
sys.path.insert(0, CTRL)
sys.argv = [SCRIPT] + sys.argv[1:]
runpy.run_path(SCRIPT, run_name='__main__')
