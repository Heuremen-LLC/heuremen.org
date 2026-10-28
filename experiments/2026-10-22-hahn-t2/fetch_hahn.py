"""fetch_hahn.py CHIP -- save raw_hahn_<chip>.npz: 'c_<pub>' as (shots, clbits) booleans, column j = classical bit j (to_bool_array, little)."""
import sys, json, numpy as np
chip = sys.argv[1]; sys.argv = [sys.argv[0], chip]
import hahn as R
rec = json.load(open(f'{R.HERE}/job_hahn_{chip}.json')); res = R.service().job(rec['job_id']).result(); out = {}
for i, pr in enumerate(res): out[f'c_{i}'] = pr.data.c.to_bool_array(order='little')
np.savez_compressed(f'{R.HERE}/raw_hahn_{chip}.npz', **out); print('saved', chip, len(res), 'pubs')
