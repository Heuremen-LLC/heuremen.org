"""fire_ladder_r2.py ARGS -- runs the UNCHANGED ladder.py source (its hash in CODE_LADDER.sha256 holds) with every saved file name given the
suffix _r2, so the repeat never overwrites run 1. Same circuits, same seed, same shots, same analysis."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); src = open(f'{HERE}/ladder.py').read()
n1 = src.count('ladder_{CHIP}.'); n2 = src.count('ladder_{chip}.'); assert n1 >= 2 and n2 >= 3, (n1, n2)
src = src.replace('ladder_{CHIP}.', 'ladder_{CHIP}_r2.').replace('ladder_{chip}.', 'ladder_{chip}_r2.')
sys.argv = ['ladder.py'] + sys.argv[1:]; sys.path.insert(0, HERE)
exec(compile(src, f'{HERE}/ladder.py', 'exec'), {'__name__': '__main__', '__file__': f'{HERE}/ladder.py'})
