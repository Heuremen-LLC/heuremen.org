"""Tabulate synth_results_full.txt (flipcheck_run2) into the numbers quoted in PLAN_FLIP.md. Heuremen LLC, Wayfinder and Bones, 22 Oct 2026."""
import re, sys, numpy as np, collections
f = sys.argv[1] if len(sys.argv) > 1 else 'synth_results_full.txt'
sc = None; rows = {}; fin = {}; dec = collections.Counter(); kg = collections.Counter(); kgf = []; idle = []
for l in open(f):
    m = re.match(r'#+ SCENARIO (\S+)', l)
    if m: sc = m.group(1); g = 'null' if sc.startswith('null') else sc[0]; continue
    m = re.match(r'\s+N=(\d) (F2h|F2c|F2|L1w|L1|L2|L3)\s+.*pooled\s+([+-][\d.]+)%?\s*\[([+-][\d.]+), ([+-][\d.]+)\]', l)
    if m: rows[(sc, m.group(2), m.group(1))] = [float(m.group(i)) for i in (3, 4, 5)]
    m = re.match(r'\s+N=(\d) (F2h|F2c|F2|L1|L2|L3)\s+FINAL -> (.*)', l)
    if m: fin[(sc, m.group(2), m.group(1))] = m.group(3).strip()[:11]
    if 'DECISION' in l:
        w = l.split('| was:')[-1]; dec[(g, ('KGFAIL+' if 'KNOWN-GOOD FAILED' in l else '') + ('EXCL' if 'EXCLUDED' in w else 'MAY' if 'MAY' in w else 'other'))] += 1
    if 'KNOWN-GOOD' in l and 'SIDE' not in l and 'DECISION' not in l:
        k = '0X' if '0X:' in l else 'idle' if 'idle N=' in l else 'L1w'; r = 'PASS' if ('PASS' in l.split('|')[0]) else 'FAIL'
        kg[(g, k, r)] += 1
        if k == 'idle': idle.append(float(re.search(r'floor ([+-][\d.]+)%', l).group(1)))
        if r == 'FAIL' and k != 'L1w': kgf.append((sc, l.strip()[:60]))
grp = {'a': ['a'], 'b': ['b'] + ['b%d' % k for k in range(1, 6)], 'null': ['null%d' % k for k in range(6)]}
for g, scs in grp.items():
    for r in ('L1', 'L2', 'L3', 'F2', 'F2h', 'F2c'):
        v = [rows[(s, r, N)] for s in scs for N in '48']; fc = collections.Counter(fin[(s, r, N)] for s in scs for N in '48')
        print('%-4s %-3s medians %+.4f..%+.4f  half-width %.4f  finals %s' % (g, r, min(x[0] for x in v), max(x[0] for x in v), np.mean([(x[2] - x[1]) / 2 for x in v]), dict(fc)))
print('decisions:', dict(dec)); print('known-good:', dict(sorted(kg.items()))); print('instrument known-good FAILs:', kgf)
idle = np.array(idle); print('idle minus unheralded 0X floor (points): n=%d mean %+.4f sd %.4f range %+.4f..%+.4f' % (len(idle), idle.mean(), idle.std(), idle.min(), idle.max()))
