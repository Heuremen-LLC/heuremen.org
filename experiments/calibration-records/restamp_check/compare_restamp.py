"""Compare T1 and T2 qubit by qubit across each restamped day pair: same value AND same parameter date = the same measurement carried
under a new record stamp; same value with a new date, or a new value, = remeasured. Reads the saved records beside this script."""
import json, gzip, glob, os
HERE = os.path.dirname(os.path.abspath(__file__))
PAIRS = [('ibm_boston','2026-09-18','2026-09-19'),('ibm_boston','2026-09-21','2026-09-22'),('ibm_boston','2026-09-23','2026-09-24'),
         ('ibm_boston','2026-09-25','2026-09-26'),('ibm_boston','2026-10-02','2026-10-03'),('ibm_marrakesh','2026-09-20','2026-09-21'),
         ('ibm_phoenix','2026-10-05','2026-10-06'),('ibm_pittsburgh','2026-09-16','2026-09-17')]
def load(chip, d):
    r = json.load(gzip.open(os.path.join(HERE, f'{chip}_{d}.json.gz'), 'rt'))
    out = {}
    for q, ps in enumerate(r['qubits']):
        for p in ps:
            if p['name'] in ('T1', 'T2'): out[(q, p['name'])] = (p['value'], p['date'])
    return r['last_update_date'], out
for chip, d0, d1 in PAIRS:
    s0, a = load(chip, d0); s1, b = load(chip, d1)
    keys = sorted(set(a) & set(b))
    same = sum(1 for k in keys if a[k] == b[k]); newdate_samevalue = sum(1 for k in keys if a[k][0] == b[k][0] and a[k][1] != b[k][1])
    changed = sum(1 for k in keys if a[k][0] != b[k][0])
    newest = max((b[k][1] for k in keys), default=None)
    print(f'{chip} {d0} -> {d1}  stamps {s0} -> {s1}  T1/T2 entries {len(keys)}: identical value and date {same}, '
          f'same value new date {newdate_samevalue}, new value {changed}; newest T1/T2 date in the later record {newest}')
