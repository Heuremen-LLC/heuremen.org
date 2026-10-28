"""fetch_iso.py CHIP -- save raw_iso_<chip>.npz: 'c_<pub>' (shots, clbits) booleans, column j = j-th listed qubit (to_bool_array, little)."""
import sys, json, numpy as np
chip = sys.argv[1]; sys.argv = [sys.argv[0], chip]
import readout_audit as R
rec = json.load(open(f'{R.HERE}/job_iso_{chip}.json')); res = R.service().job(rec['job_id']).result(); out = {}
for i, pr in enumerate(res): out[f'c_{i}'] = pr.data.c.to_bool_array(order='little')
np.savez_compressed(f'{R.HERE}/raw_iso_{chip}.npz', **out); print('saved', chip, len(res), 'pubs')
