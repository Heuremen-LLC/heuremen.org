"""fetch_hahn_disco.py CHIP -- fetch_hahn.py for jobs fired through fire_hahn_disco.py (Darkwave Disco instance)."""
import sys, json, numpy as np
chip = sys.argv[1]; sys.argv = [sys.argv[0], chip]
import hahn as R
R.INSTANCE = 'Darkwave Disco'
rec = json.load(open(f'{R.HERE}/job_hahn_{chip}.json')); res = R.service().job(rec['job_id']).result(); out = {}
for i, pr in enumerate(res): out[f'c_{i}'] = pr.data.c.to_bool_array(order='little')
np.savez_compressed(f'{R.HERE}/raw_hahn_{chip}.npz', **out); print('saved', chip, len(res), 'pubs')
