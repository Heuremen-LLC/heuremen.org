"""fetch_threebox.py CHIP -- pull the job's result and save raw_<chip>.npz keyed '<var>|<box>|<eps>|<copy>' (packed bits: bit0 s0, bit1 s1, bit2 probe)."""
import sys, json, os, numpy as np
import run_threebox as R
HERE = R.HERE
CHIP = sys.argv[1]; meta = json.load(open(f'{HERE}/job_{CHIP}.json')); R.KEYFILE, R.INSTANCE = meta['keyfile'], meta['instance']; job = R.service().job(meta['job']); res = job.result(); out = {}
for (kk, cp), r in zip(meta['items'], res):
    out[f'{kk[0]}|{kk[1]}|{kk[2]:.2f}|{cp}'] = np.asarray(r.data.c.array[:, 0]).astype(np.uint8)
np.savez_compressed(f'{HERE}/raw_{CHIP}.npz', **out); print('saved', len(out), 'circuits for', CHIP, meta['job'])
