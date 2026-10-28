"""fire_hahn_key2.py CHIP [--fire] -- runs the unchanged hahn.py source with the key swapped to the open instance of
.ibm_quantum_key_2, per ADDENDUM_HAHN_FRESH_KEY2 (written before any fresh-record job). hahn.py itself is not edited,
so CODE_HAHN.sha256 still holds."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
src = open(f'{HERE}/hahn.py').read()
old = "KEYFILE, INSTANCE = '.ibm_quantum_key_3', 'open-instance'"
assert src.count(old) == 1
src = src.replace(old, "KEYFILE, INSTANCE = '.ibm_quantum_key_2', 'open-instance'")
sys.argv = ['hahn.py'] + sys.argv[1:]
g = {'__name__': '__main__', '__file__': f'{HERE}/hahn.py'}
exec(compile(src, f'{HERE}/hahn.py', 'exec'), g)
