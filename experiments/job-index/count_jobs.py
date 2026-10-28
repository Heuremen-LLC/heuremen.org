"""Counts behind the paper's run table. python3 count_jobs.py (reads the two JSON files beside it)."""
import json, collections, datetime as dt
ibm = json.load(open('ibm_jobs_2026-09-24_to_10-22.json'))
print('IBM jobs', len(ibm), dict(collections.Counter(j['status'] for j in ibm)))
day = collections.defaultdict(collections.Counter)
for j in ibm:
    day[j['created_et'][:10]][j['backend']] += 1
for d in sorted(day):
    print(' ', d, 'US Eastern', sum(day[d].values()), dict(sorted(day[d].items())))
rig = json.load(open('rigetti_tasks_2026-09-25_to_10-05.json'))
et = lambda s: (dt.datetime.fromisoformat(s) - dt.timedelta(hours=4)).date().isoformat()
print('Rigetti tasks', len(rig), dict(collections.Counter(t['status'] for t in rig)))
rd = collections.defaultdict(collections.Counter)
for t in rig:
    rd[et(t['created'])][t['status']] += 1
for d in sorted(rd):
    print(' ', d, 'US Eastern', sum(rd[d].values()), dict(rd[d]))
used = sum(1 for t in rig if t['created'][:10] == '2026-10-04' and t['created'][11:13] in ('16', '17'))
print('Rigetti 4 October, 16 to 17 UTC (used in the paper):', used)
