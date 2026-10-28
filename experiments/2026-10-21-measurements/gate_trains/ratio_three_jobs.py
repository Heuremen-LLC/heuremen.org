"""Start-|1> over start-|0> excess against one exponential at 5, 10, 20 us on 4 October, every chip, with pittsburgh's three jobs."""
import os
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'verify_xtrain.py')).read().split("def main():")[0]
exec(src)
import re
runs = {'ibm_marrakesh': (301, 302, 303), 'ibm_kingston': (301, 302, 303), 'ibm_pittsburgh': (361, 362, 363)}
allr = []
for chip, rr in runs.items():
    rs = sweep_chip(chip, rr)
    for t in (5, 10, 20):
        m0, _ = sw(rs, 'OFF', t, 'ONE'); m1, _ = sw(rs, 'ON', t, 'ONE'); r = -m1 / m0; allr.append(r)
        print(f'{chip:<15} {t:>3} us  start0 {m0:+.4f}  start1 {m1:+.4f}  ratio {100*r:.0f} %')
print(f'range {100*min(allr):.0f} to {100*max(allr):.0f} %')
