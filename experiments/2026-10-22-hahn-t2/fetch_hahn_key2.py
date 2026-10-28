"""fetch_hahn_key2.py CHIP -- fetch_hahn.py with hahn's KEYFILE swapped to .ibm_quantum_key_2 (open instance), per ADDENDUM_HAHN_FRESH_KEY2."""
import sys, json, numpy as np, importlib.util
chip = sys.argv[1]; sys.argv = [sys.argv[0], chip]
HERE = '/Users/heuremen/bones-forge/arxiv_one_readout/hahn_t2'
src = open(f'{HERE}/hahn.py').read(); old = "KEYFILE, INSTANCE = '.ibm_quantum_key_3', 'open-instance'"; assert src.count(old) == 1
spec = importlib.util.spec_from_loader('hahn', loader=None); R = importlib.util.module_from_spec(spec); R.__file__ = f'{HERE}/hahn.py'
exec(compile(src.replace(old, "KEYFILE, INSTANCE = '.ibm_quantum_key_2', 'open-instance'"), f'{HERE}/hahn.py', 'exec'), R.__dict__)
rec = json.load(open(f'{HERE}/job_hahn_{chip}.json')); res = R.service().job(rec['job_id']).result(); out = {}
for i, pr in enumerate(res): out[f'c_{i}'] = pr.data.c.to_bool_array(order='little')
np.savez_compressed(f'{HERE}/raw_hahn_{chip}.npz', **out); print('saved', chip, len(res), 'pubs')
