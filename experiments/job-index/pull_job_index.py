"""Re-pull of every IBM job on our credentials, 23 Oct 2026 (the 22 Oct file was pulled at 16:28 EDT and missed that evening's jobs).
Every key file, every instance it reaches; job ID, processor, creation time (UTC and US Eastern), status. Writes
ibm_jobs_2026-09-24_to_10-22.json (US Eastern creation day 24 Sep to 22 Oct, as the paper's period) and ibm_jobs_2026-10-23.json (later)."""
import json, os, datetime as dt
from zoneinfo import ZoneInfo
from qiskit_ibm_runtime import QiskitRuntimeService
H = os.path.expanduser('~/bones-forge'); ET = ZoneInfo('America/New_York'); seen = {}
SINCE = dt.datetime(2026, 9, 24, 4, tzinfo=dt.timezone.utc)
for f in ['.ibm_quantum_key', '.ibm_quantum_key_2', '.ibm_quantum_key_3', '.ibm_quantum_key_4', '.ibm_quantum_key_parallax']:
    p = os.path.join(H, f)
    if not os.path.exists(p): print(f, 'absent'); continue
    tok = open(p).read().strip()
    try:
        base = QiskitRuntimeService(channel='ibm_quantum_platform', token=tok)
        for inst in base.instances():
            s = QiskitRuntimeService(channel='ibm_quantum_platform', token=tok, instance=inst['crn']); n = 0
            for j in s.jobs(limit=None, created_after=SINCE):
                if j.job_id() in seen: continue
                c = j.creation_date.astimezone(dt.timezone.utc); n += 1
                seen[j.job_id()] = dict(job_id=j.job_id(), backend=j.backend().name if j.backend() else None, created_utc=c.isoformat(),
                                        created_et=c.astimezone(ET).isoformat(), status=str(j.status()), pubs=None, tags=list(j.tags or []) or None)
            print(f, inst.get('name'), 'new', n, flush=True)
    except Exception as e: print(f, 'no:', type(e).__name__)
rows = sorted(seen.values(), key=lambda r: r['created_utc'])
a = [r for r in rows if r['created_et'][:10] <= '2026-10-22']; b = [r for r in rows if r['created_et'][:10] > '2026-10-22']
json.dump(a, open('ibm_jobs_2026-09-24_to_10-22.json', 'w'), indent=0); json.dump(b, open('ibm_jobs_2026-10-23.json', 'w'), indent=0)
print('total', len(rows), 'to 22 Oct', len(a), 'from 23 Oct', len(b))
