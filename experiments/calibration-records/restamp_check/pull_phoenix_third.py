"""Third phoenix record of 6 Oct (12:25 EDT), named in the daily-medians row as identical to 5 Oct 09:44 and 6 Oct 10:45: pulled from IBM's
properties history (no job) and compared qubit by qubit with the 6 Oct 10:45 record already in this folder."""
import json, os, sys, gzip, datetime as D
sys.path.insert(0, os.path.expanduser('~/bones-forge/parallax/out/qubit_fifth_run'))
import live
HERE = os.path.dirname(os.path.abspath(__file__))
svc, _ = live.open_service('RESTAMP CHECK: phoenix 6 Oct 12:25 record; no job')
p = svc.backend('ibm_phoenix').properties(datetime=D.datetime.fromisoformat('2026-10-06T16:40:00+00:00'))
f = os.path.join(HERE, 'ibm_phoenix_2026-10-06_1225.json.gz'); json.dump(p.to_dict(), gzip.open(f, 'wt'), default=str)
print('saved', p.last_update_date)
