"""shot_grid.py -- which shot grid IBM's published readout probabilities sit on, per saved record (23 Oct 2026, written after row 9 of
QIA-2026-001 was drafted). A value that is a multiple of 1/4096 comes from 4096 shots per state; all values on 1/1024 means 1024 shots."""
import json, gzip, glob
for chip in ('ibm_fez', 'ibm_kingston', 'ibm_marrakesh'):
    for f in sorted(glob.glob(f'**/*{chip}*.json*', recursive=True)):
        try: d = json.load(gzip.open(f) if f.endswith('.gz') else open(f))
        except Exception: continue
        if not isinstance(d, dict) or 'qubits' not in d: continue
        v = [x['value'] for q in d['qubits'] for x in q if x['name'] in ('prob_meas1_prep0', 'prob_meas0_prep1') and x['value'] not in (None, 0, 1)]
        if not v: continue
        g = lambda n: sum(abs(x * n - round(x * n)) < 1e-6 for x in v) / len(v)
        print(f'{chip} record {str(d.get("last_update_date"))[:16]}: on 1/1024 grid {g(1024):.2f}, on 1/4096 grid {g(4096):.2f} ({len(v)} values) -> {"1024 shots" if g(1024) > 0.99 else "4096 shots"}')
