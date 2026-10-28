"""Fire direct.py for ibm_phoenix on Darkwave Disco, as PLAN_DIRECT.md says ("ibm_phoenix on Darkwave Disco when its usage cap allows;
phoenix is not on the open instance"). direct.py hard-codes the open instance; this wrapper runs its unchanged source with only the
key file and instance name swapped, so direct.py's hash stays as committed. Usage: python3 fire_phoenix_disco.py [--fire]"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(HERE, 'direct.py')).read()
old = "KEYFILE, INSTANCE = '.ibm_quantum_key', 'open-instance'"
assert src.count(old) == 1
src = src.replace(old, "KEYFILE, INSTANCE = '.ibm_quantum_key_3', 'Darkwave Disco'")
sys.argv = [os.path.join(HERE, 'direct.py'), 'ibm_phoenix'] + (['--fire'] if '--fire' in sys.argv else [])
g = {'__name__': '__main__', '__file__': os.path.join(HERE, 'direct.py')}
exec(compile(src, os.path.join(HERE, 'direct.py'), 'exec'), g)
