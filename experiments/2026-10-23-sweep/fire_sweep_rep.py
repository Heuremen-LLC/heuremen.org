"""fire_sweep_rep.py TAG ARGS -- runs the UNCHANGED sweep.py source (hash in CODE_SWEEP.sha256) with every saved file name suffixed _TAG,
so a repeat never overwrites run 1. Same circuits, same seed, same shots, same analysis."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); tag = sys.argv[1]; src = open(f'{HERE}/sweep.py').read()
n = sum(src.count(f'_sweep_{{{v}}}.') for v in ('CHIP', 'chip')); assert n >= 6, n
for v in ('CHIP', 'chip'): src = src.replace(f'_sweep_{{{v}}}.', f'_sweep_{{{v}}}_{tag}.')
sys.argv = ['sweep.py'] + sys.argv[2:]; sys.path.insert(0, HERE)
exec(compile(src, f'{HERE}/sweep.py', 'exec'), {'__name__': '__main__', '__file__': f'{HERE}/sweep.py'})
