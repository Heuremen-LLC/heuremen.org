"""hahn_table.py -- every Hahn T2 run side by side (defined after the data, 23 Oct 2026 morning; reporting only, no new fit).
Per run: job, when IBM ran it (job metrics), the record saved at submission, the median age of the stated T2 values at run time,
and the single-run reading under both fits of hahn_pool.py (FIXED 1/2, FREE asymptote). Then the clock (ADDENDUM_HAHN_CLOCK_FLOOR.md):
decline from the earliest to the latest fez point against twice the largest between-job shift at near-equal age, under both fits."""
import sys, json, subprocess, statistics as st, datetime as dt, re
sys.path.insert(0, '.'); import hahn as R
RUNS = [('ibm_kingston','run 1','kingston_run1'),('ibm_kingston','run 2','kingston_run2'),('ibm_marrakesh','run 1','marrakesh_run1'),
        ('ibm_marrakesh','run 2','marrakesh_run2'),('ibm_fez','run 1','fez_run1'),('ibm_fez','run 2','fez_run2'),('ibm_fez','clock 1','fez_clock1'),
        ('ibm_fez','clock 2','fez_clock2'),('ibm_pittsburgh','run 1','pittsburgh'),('ibm_phoenix','run 1','phoenix'),
        ('ibm_kingston','fresh record','kingston_fresh'),('ibm_boston','run 1','boston')]
PAT = re.compile(r'median ([\d.]+) \[([\d.]+), ([\d.]+)\] over (\d+) qubits -> (\w+)\s+\(below stated (\d+)\)')
val = {}
for chip, run, tag in RUNS:
    job = json.load(open(f'job_hahn_ibm_{tag}.json')); props = json.load(open(f'props_hahn_ibm_{tag}.json'))
    R.INSTANCE = 'Darkwave Disco' if job['instance'] != 'open-instance' else 'open-instance'
    R.KEYFILE = '.ibm_quantum_key_2' if tag.endswith('fresh') else '.ibm_quantum_key_3'   # 23 Oct: fresh-record runs fired on key 2 (ADDENDUM_HAHN_FRESH_KEY2)
    ts = R.service().job(job['job_id']).metrics().get('timestamps', {}); ran = dt.datetime.fromisoformat(ts['running'].replace('Z', '+00:00'))
    ages = [(ran - dt.datetime.fromisoformat(x['date'])).total_seconds() / 3600 for q in props['qubits'] for x in q if x['name'] == 'T2']
    pool = subprocess.run([sys.executable, 'hahn_pool.py', chip, f'raw_hahn_ibm_{tag}.npz', f'job_hahn_ibm_{tag}.json'], capture_output=True, text=True).stdout
    fx, fr = [PAT.search(l).groups() for l in pool.splitlines() if 'FIXED' in l or 'FREE' in l]
    val[tag] = (st.median(ages), float(fx[0]), float(fr[0]))
    print(f"{chip} {run} job {job['job_id']} ran {ran:%Y-%m-%d %H:%M} UTC; record {props['last_update_date']}; stated T2 median age at run {st.median(ages):.1f} h; queued {(ran - dt.datetime.fromisoformat(job['submitted_utc'].replace('Z', '+00:00'))).total_seconds() / 3600:.1f} h")
    print(f"  FIXED {fx[0]} [{fx[1]}, {fx[2]}] -> {fx[4]} ({fx[5]} of {fx[3]} below stated) | FREE {fr[0]} [{fr[1]}, {fr[2]}] -> {fr[4]} ({fr[5]} of {fr[3]} below stated)")
T2 = lambda t: {(q, str(x['date']), x['value']) for q, ps in enumerate(json.load(open(f'props_hahn_ibm_{t}.json'))['qubits']) for x in ps if x['name'] == 'T2'}
same = set.intersection(*[T2(t) for t in ('fez_run1', 'fez_run2', 'fez_clock1', 'fez_clock2')])
print(f'\nCLOCK (ibm_fez, stated T2 values and dates identical in all four records, {len(same)} of {len(T2("fez_run1"))}):')
pts = ['fez_run1', 'fez_run2', 'fez_clock1', 'fez_clock2']
for k, name in ((1, 'FIXED'), (2, 'FREE')):
    first, last = val[pts[0]][k], val[pts[-1]][k]
    spread = max(abs(val['fez_run1'][k] - val['fez_run2'][k]), abs(val['kingston_run1'][k] - val['kingston_run2'][k]), abs(val['marrakesh_run1'][k] - val['marrakesh_run2'][k]))
    out = 'BELOW 0 (decays)' if first - last > 2 * spread else 'CONSISTENT'
    print(f"  {name}: {' -> '.join(f'{val[p][k]:.3f}' for p in pts)} at {' / '.join(f'{val[p][0]:.1f}' for p in pts)} h; decline first to last {first - last:+.3f}; floor 2 x {spread:.3f} = {2 * spread:.3f} -> {out}")
print('  Points at 8 and 16 h were not run.')
d = {t: abs(v[1] - v[2]) for t, v in val.items()}; t = max(d, key=d.get)
print(f'FIT MODEL: largest |planned - free| over the twelve runs {d[t]:.3f} ({t})')
