import json, sys
R = json.load(open('results_oct07.json'))
def cell(p):
    if p is None: return 'copies shifted -> CANNOT TELL'
    m, lo, hi, n, o = p
    if m != m: return 'n=%d -> CANNOT TELL' % n
    return '%+.3f [%+.3f,%+.3f] n=%d %s' % (m, lo, hi, n, o.replace(' ', '_'))
rows = {}
for c, f, l, k, p in R: rows.setdefault((c, f, l), {})[k] = p
print('SUMMARY TABLE (paper convention: kick negative = same direction as 21-22 Oct; deg per readout; bins = couplings to nearest looked bit)')
for (c, f, l), d in rows.items():
    if 'kick all' not in d: continue
    print('%s | %s | %s' % (c.replace('ibm_', ''), f, l))
    for k in ('kick all', 'kick d<=2', 'kick d=3-4', 'kick d>=5', 'kick d>=3', 'kickN4', 'coh all', 'coh d<=2', 'coh d>=3', 'cohN4'):
        if k in d: print('    %-11s %s' % (k, cell(d[k])))
