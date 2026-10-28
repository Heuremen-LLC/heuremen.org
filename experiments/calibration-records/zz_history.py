"""zz of five live edges across every IBM calibration record published on this site (this folder's September and October records and
the records saved beside each October run). IBM reports zz in GHz; printed here in kHz. Standard library only.
  python3 zz_history.py > results_zz_history.text"""
import gzip, json, os
EXP = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EDGES = {'ibm_fez': ['zz_1314', 'zz_3953'], 'ibm_phoenix': ['zz_100101', 'zz_7080', 'zz_8696']}
seen = {}
for dp, dn, fn in os.walk(EXP):
    for f in fn:
        if not (f.endswith('.json') or f.endswith('.json.gz')): continue
        p = os.path.join(dp, f)
        if os.path.getsize(p) > 30e6: continue
        try: d = json.load(gzip.open(p, 'rt') if f.endswith('.gz') else open(p))
        except Exception: continue
        if not (isinstance(d, dict) and d.get('backend_name') in EDGES and 'general' in d): continue
        g = {x['name']: x for x in d['general']}
        for k in EDGES[d['backend_name']]:
            if k in g: seen.setdefault((d['backend_name'], k), {})[d['last_update_date'][:16]] = (g[k]['value'] * 1e6, os.path.relpath(p, EXP))
for (chip, k), recs in sorted(seen.items()):
    print(f'{chip} {k}')
    for stamp in sorted(recs): print(f'   record {stamp}  zz {recs[stamp][0]:6.2f} kHz   ({recs[stamp][1]})')
