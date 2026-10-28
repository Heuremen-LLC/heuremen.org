"""fetch_direct_disco.py CHIP -- the raw-saving part of fetch_direct.py for a job fired through fire_phoenix_disco.py (Darkwave Disco).
The calibration record was saved AT FIRING (props_direct_<chip>.json), so the after-the-fact history pull is not repeated."""
import sys, json, numpy as np
chip = sys.argv[1]; sys.argv = [sys.argv[0], chip]
import direct as R
R.KEYFILE, R.INSTANCE = '.ibm_quantum_key_3', 'Darkwave Disco'
rec = json.load(open(f'{R.HERE}/job_direct_{chip}.json')); job = R.service().job(rec['job_id']); print('status', job.status())
res = job.result(); out = {}
for i, pr in enumerate(res):
    for name in pr.data: out['%s_%d' % (name, i)] = getattr(pr.data, name).array
np.savez_compressed(f'{R.HERE}/raw_direct_{chip}.npz', **out); print('saved raw', len(res), 'pubs')
