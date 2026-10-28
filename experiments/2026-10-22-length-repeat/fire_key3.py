"""fire_key3.py CHIP [--shots N] [--fire] -- runs the UNCHANGED run_len.py source with KEYFILE swapped to .ibm_quantum_key_3 (open instance)."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); src = open(os.path.join(HERE, 'run_len.py')).read()
old = "KEYFILE, INSTANCE = '.ibm_quantum_key', 'open-instance'"; assert src.count(old) == 1
src = src.replace(old, "KEYFILE, INSTANCE = '.ibm_quantum_key_3', 'open-instance'")
sys.argv = [os.path.join(HERE, 'run_len.py')] + sys.argv[1:]
exec(compile(src, os.path.join(HERE, 'run_len.py'), 'exec'), {'__name__': '__main__', '__file__': os.path.join(HERE, 'run_len.py')})
