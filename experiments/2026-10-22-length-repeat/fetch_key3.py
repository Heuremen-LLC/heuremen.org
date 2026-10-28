"""fetch_key3.py CHIP -- fetch_len.py for a job fired through fire_key3.py (KEYFILE .ibm_quantum_key_3)."""
import sys, json, numpy as np
chip = sys.argv[1]; sys.argv = [sys.argv[0], chip]
import run_len as R
R.KEYFILE = '.ibm_quantum_key_3'
exec(open(R.HERE + '/fetch_len.py').read().split('import run_len as R', 1)[1])
