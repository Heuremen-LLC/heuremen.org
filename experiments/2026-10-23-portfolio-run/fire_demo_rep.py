"""fire_demo_rep.py TAG ARGS -- runs the UNCHANGED demo.py source (hash in CODE_DEMO.sha256) with saved file names suffixed _TAG."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); tag = sys.argv[1]; src = open(f'{HERE}/demo.py').read()
n = src.count('_demo_{CHIP}.'); assert n >= 5, n
src = src.replace('_demo_{CHIP}.', f'_demo_{{CHIP}}_{tag}.')
sys.argv = ['demo.py'] + sys.argv[2:]; sys.path.insert(0, HERE)
exec(compile(src, f'{HERE}/demo.py', 'exec'), {'__name__': '__main__', '__file__': f'{HERE}/demo.py'})
