"""fetch_len.py CHIP -- save raw_<chip>.npz: '<register>_<pub index>' packed bit arrays, as the control run."""
import sys, json, numpy as np
sys.argv = [sys.argv[0], sys.argv[1]]
import run_len as R
rec = json.load(open(f'{R.HERE}/job_{sys.argv[1]}.json')); res = R.service().job(rec['job_id']).result(); out = {}
for i, pr in enumerate(res):
    for name in pr.data: out['%s_%d' % (name, i)] = getattr(pr.data, name).array
np.savez_compressed(f'{R.HERE}/raw_{sys.argv[1]}.npz', **out); print('saved', sys.argv[1], len(res), 'pubs')
