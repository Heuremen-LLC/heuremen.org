"""fire_hahn_disco.py CHIP [--fire] -- runs the unchanged hahn.py source with the instance swapped to Darkwave Disco (paid),
for chips the open instance does not reach. hahn.py itself is not edited, so CODE_HAHN.sha256 still holds. Also prints the median
age of IBM's stated T2 at the moment of firing."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
src = open(f'{HERE}/hahn.py').read()
old = "KEYFILE, INSTANCE = '.ibm_quantum_key_3', 'open-instance'"
assert src.count(old) == 1
src = src.replace(old, "KEYFILE, INSTANCE = '.ibm_quantum_key_3', 'Darkwave Disco'")
sys.argv = ['hahn.py'] + sys.argv[1:]
g = {'__name__': '__main__', '__file__': f'{HERE}/hahn.py'}
exec(compile(src, f'{HERE}/hahn.py', 'exec'), g)
