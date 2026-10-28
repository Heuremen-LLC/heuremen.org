"""Same comparison as compare_restamp.py (T1 and T2 per qubit: value AND parameter date) on the record pairs of 22 and 23 Oct 2026
saved by pull_oct23.py. Every pair is printed, whatever it shows."""
import json, gzip, os
HERE = os.path.dirname(os.path.abspath(__file__))
PAIRS = [('ibm_fez', '2026-10-22_1817', '2026-10-22_1955'), ('ibm_fez', '2026-10-22_1955', '2026-10-22_2149'), ('ibm_fez', '2026-10-22_2149', '2026-10-23_0438'),
         ('ibm_kingston', '2026-10-22_1845', '2026-10-23_0455'), ('ibm_marrakesh', '2026-10-22_1821', '2026-10-23_0415')]
def load(chip, d):
    r = json.load(gzip.open(os.path.join(HERE, f'{chip}_{d}.json.gz'), 'rt')); out = {}
    for q, ps in enumerate(r['qubits']):
        for p in ps:
            if p['name'] in ('T1', 'T2'): out[(q, p['name'])] = (p['value'], str(p['date']))
    return r['last_update_date'], out
for chip, d0, d1 in PAIRS:
    s0, a = load(chip, d0); s1, b = load(chip, d1); keys = sorted(set(a) & set(b))
    same = sum(a[k] == b[k] for k in keys); nd = sum(a[k][0] == b[k][0] and a[k][1] != b[k][1] for k in keys); nv = sum(a[k][0] != b[k][0] for k in keys)
    by = {n: sum(a[k] == b[k] for k in keys if k[1] == n) for n in ('T1', 'T2')}; tot = {n: sum(k[1] == n for k in keys) for n in ('T1', 'T2')}
    print(f'{chip} stamps {s0} -> {s1}  T1/T2 entries {len(keys)}: identical value and date {same} (T1 {by["T1"]} of {tot["T1"]}, T2 {by["T2"]} of {tot["T2"]}), '
          f'same value new date {nd}, new value {nv}; newest T1/T2 date in the later record {max(b[k][1] for k in keys)}')
import re
old = open(os.path.join(HERE, 'results_restamp.text')).read()
o_same = sum(int(x) for x in re.findall(r'identical value and date (\d+)', old)); o_all = sum(int(x) for x in re.findall(r'T1/T2 entries (\d+)', old))
n_same = n_all = 0
for chip, d0, d1 in PAIRS:
    _, a = load(chip, d0); _, b = load(chip, d1); keys = set(a) & set(b); n_all += len(keys); n_same += sum(a[k] == b[k] for k in keys)
print(f'TOTAL: earlier pairs {o_same} of {o_all} (results_restamp.text), these pairs {n_same} of {n_all}; all {o_same + n_same:,} of {o_all + n_all:,} entries identical in value and date across {len(re.findall("stamps", old)) + len(PAIRS)} record pairs')
print('\nEVERY CONSECUTIVE PAIR IN THIS FOLDER (by record stamp, per chip), counted or not:')
import collections, glob, datetime as dt
recs = collections.defaultdict(dict)
for f in sorted(glob.glob(os.path.join(HERE, 'ibm_*.json.gz'))):
    n = os.path.basename(f)[:-8]; chip = '_'.join(n.split('_')[:2]); st, ent = load(chip, n[len(chip) + 1:]); recs[chip][str(st)] = ent
for chip in sorted(recs):
    ss = sorted(recs[chip], key=lambda x: dt.datetime.fromisoformat(x))
    for s0, s1 in zip(ss, ss[1:]):
        a, b = recs[chip][s0], recs[chip][s1]; keys = set(a) & set(b); gap = (dt.datetime.fromisoformat(s1) - dt.datetime.fromisoformat(s0)).total_seconds() / 3600
        print(f'  {chip} {s0} -> {s1}  gap {gap:.1f} h  identical {sum(a[k] == b[k] for k in keys)} of {len(keys)}')
print('\nWHAT ELSE CHANGED UNDER THE SAME STAMPS (gate_error per gate and qubits, readout_error per qubit):')
def other(chip, d):
    r = json.load(gzip.open(os.path.join(HERE, f'{chip}_{d}.json.gz'), 'rt'))
    g = {(x['gate'], tuple(x['qubits'])): p['value'] for x in r.get('gates', []) for p in x['parameters'] if p['name'] == 'gate_error'}
    ro = {i: p['value'] for i, q in enumerate(r['qubits']) for p in q if p['name'] == 'readout_error'}
    return g, ro
for chip, d0, d1 in PAIRS + [('ibm_phoenix', '2026-10-05', '2026-10-06'), ('ibm_phoenix', '2026-10-06', '2026-10-06_1225')]:
    (ga, ra), (gb, rb) = other(chip, d0), other(chip, d1); ch = collections.defaultdict(set)
    for k in ga:
        if k in gb and ga[k] != gb[k]: ch[k[0]].add(k[1])
    print(f'  {chip} {d0} -> {d1}: gate_error changed ' + (', '.join(f'{g} on {len(v)}' for g, v in sorted(ch.items())) or 'none') + f'; readout_error changed on {sum(ra[i] != rb.get(i) for i in ra)} qubits')
