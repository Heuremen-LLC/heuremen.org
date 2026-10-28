"""fetch_direct.py CHIP -- read the submitted job (no new job): save raw_direct_<chip>.npz ('<register>_<pub index>' packed bit arrays,
as the other runs) and props_direct_<chip>.json, IBM's calibration in force at the job's submission time, pulled from IBM's properties
history after the fact (it was not saved at firing)."""
import sys, os, json, datetime as D
import numpy as np
chip = sys.argv[1]; sys.argv = [sys.argv[0], chip]
import direct as R
rec = json.load(open(f'{R.HERE}/job_direct_{chip}.json')); svc = R.service(); job = svc.job(rec['job_id'])
print('status', job.status())
res = job.result(); out = {}
for i, pr in enumerate(res):
    for name in pr.data: out['%s_%d' % (name, i)] = getattr(pr.data, name).array
np.savez_compressed(f'{R.HERE}/raw_direct_{chip}.npz', **out); print('saved raw', len(res), 'pubs')
t = D.datetime.fromisoformat(rec['submitted_utc'].replace('Z', '+00:00'))
p = svc.backend(chip).properties(datetime=t)
json.dump(p.to_dict(), open(f'{R.HERE}/props_direct_{chip}.json', 'w'), default=str); print('saved props, record stamp', p.last_update_date)
