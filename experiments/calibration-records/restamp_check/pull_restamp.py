"""Per-qubit check of the 8 day pairs whose T1/T2 medians were identical under a new record stamp (fleet_records, 15 Sep to 6 Oct 2026).
Pulls IBM's published properties for both days (same convention as the median pull: properties(datetime=<day>T12:00:00)), saves each record
as returned, and compares T1 and T2 qubit by qubit: value and the parameter's own date. No job, no QPU time."""
import json, os, sys, gzip, datetime as D
sys.path.insert(0, os.path.expanduser('~/bones-forge/parallax/out/qubit_fifth_run'))
import live
PAIRS = [('ibm_boston','2026-09-18','2026-09-19'),('ibm_boston','2026-09-21','2026-09-22'),('ibm_boston','2026-09-23','2026-09-24'),
         ('ibm_boston','2026-09-25','2026-09-26'),('ibm_boston','2026-10-02','2026-10-03'),('ibm_marrakesh','2026-09-20','2026-09-21'),
         ('ibm_phoenix','2026-10-05','2026-10-06'),('ibm_pittsburgh','2026-09-16','2026-09-17')]
HERE = os.path.dirname(os.path.abspath(__file__))
svc, _ = live.open_service('RESTAMP CHECK: published properties for 8 day pairs, per qubit; no job')
for chip, d0, d1 in PAIRS:
    be = svc.backend(chip)
    for d in (d0, d1):
        f = os.path.join(HERE, f'{chip}_{d}.json.gz')
        if os.path.exists(f): continue
        p = be.properties(datetime=D.datetime.fromisoformat(d + 'T12:00:00'))
        json.dump(p.to_dict(), gzip.open(f, 'wt'), default=str)
        print('saved', chip, d, p.last_update_date, flush=True)
