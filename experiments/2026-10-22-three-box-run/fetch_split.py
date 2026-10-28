"""fetch_split.py CHIP -- save raw_split_<chip>.npz: 'c_<pub>' as per-shot integers, bit k = classical bit k (bit0 s0, bit1 s1, bit2 probe),
the layout probe_split.xp and tb_common.xp_from_bits read."""
import sys, json, numpy as np
chip = sys.argv[1]; sys.argv = [sys.argv[0], chip]
import probe_split as R
rec = json.load(open(f'{R.HERE}/job_split_{chip}.json')); res = R.service().job(rec['job_id']).result(); out = {}
for i, pr in enumerate(res): out[f'c_{i}'] = np.array([int(s, 2) for s in pr.data.c.get_bitstrings()], dtype=np.int64)
np.savez_compressed(f'{R.HERE}/raw_split_{chip}.npz', **out); print('saved', chip, len(res), 'pubs')
