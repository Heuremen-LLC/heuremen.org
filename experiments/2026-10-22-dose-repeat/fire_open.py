"""fire_open.py [--fire] -- runs the UNCHANGED committed control-run run.py source on ibm_fez with --arms S,C, the instance name swapped
from 'Darkwave Disco' to 'open-instance' (free key), and HERE set to this folder so the 21 Oct published files are not overwritten."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(HERE, 'run.py')).read(); old = "if i.get('name') == 'Darkwave Disco'"
assert src.count(old) == 1; src = src.replace(old, "if i.get('name') == 'open-instance'")
sys.argv = [os.path.join(HERE, 'run.py'), '--chip', 'ibm_fez', '--arms', 'S,C'] + (['--fire'] if '--fire' in sys.argv else [])
exec(compile(src, os.path.join(HERE, 'run.py'), 'exec'), {'__name__': '__main__', '__file__': os.path.join(HERE, 'run.py')})
