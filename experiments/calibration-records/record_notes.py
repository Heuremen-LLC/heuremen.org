"""Record notes of QIA-2026-001 (rows 2, 3, 4, 6, 7), recomputed from IBM calibration records published on this site:
the September records in this folder and the records saved beside each October run (props_*.json, *backend_properties*.json).

  python3 record_notes.py > results_record_notes.text        (run from heuremen.org/experiments/calibration-records)

For every record: how many qubits carry each parameter, how old each value is against the record's own stamp (over 24 h, over 7 days,
oldest); and, between consecutive records of a processor, the per-qubit median relative change of readout_error, with the hours between
the two records' readout_error dates. Standard library only."""
import gzip, json, os, statistics as S
from datetime import datetime

P = datetime.fromisoformat
HERE = os.path.dirname(os.path.abspath(__file__))
EXP = os.path.dirname(HERE)


def records():
    seen = {}
    for dp, dn, fn in os.walk(EXP):
        if 'synth' in dp or 'restamp_check' in dp:
            continue
        for f in fn:
            if not (f.endswith('.json') or f.endswith('.json.gz')):
                continue
            p = os.path.join(dp, f)
            if os.path.getsize(p) > 30e6:
                continue
            try:
                d = json.load(gzip.open(p, 'rt') if f.endswith('.gz') else open(p))
            except Exception:
                continue
            if not (isinstance(d, dict) and 'qubits' in d and 'backend_name' in d and 'last_update_date' in d):
                continue
            key = (d['backend_name'], d['last_update_date'])
            seen.setdefault(key, (d, os.path.relpath(p, EXP)))
    return seen


def ages(d, par):
    st = P(d['last_update_date']); out = []; absent = 0
    for q in d['qubits']:
        v = [x for x in q if x['name'] == par and x.get('value') is not None]
        if v:
            out.append((st - P(v[0]['date'])).total_seconds() / 3600)
        else:
            absent += 1
    return out, absent


def main():
    recs = records()
    print('PER RECORD: parameter: present / absent; values over 24 h and over 7 days old against the record stamp; oldest (h)')
    for (chip, stamp), (d, src) in sorted(recs.items(), key=lambda kv: (kv[0][0], P(kv[0][1]))):
        print(f'{chip} record {stamp}  ({src})')
        for par in ('T1', 'T2', 'init_error', 'readout_error'):
            a, ab = ages(d, par)
            if a:
                print(f'   {par:<14} {len(a):>3} / {ab:<3} >24h {sum(x > 24 for x in a):>3}  >7d {sum(x > 168 for x in a):>3}  oldest {max(a):8.1f}')
            else:
                print(f'   {par:<14}   0 / {ab}')
    print('\nREADOUT ERROR BETWEEN CONSECUTIVE RECORDS (by readout_error date): per-qubit median |new/old - 1|')
    byc = {}
    for (chip, stamp), (d, src) in recs.items():
        ro = []
        for q in d['qubits']:
            v = {x['name']: x for x in q}
            if 'readout_error' not in v:
                ro = None; break
            ro.append((v['readout_error']['value'], v['readout_error']['date']))
        if ro:
            byc.setdefault(chip, {})[ro[0][1]] = [r[0] for r in ro]
    for chip in sorted(byc):
        ks = sorted(byc[chip], key=P)
        for a, b in zip(ks, ks[1:]):
            A, B = byc[chip][a], byc[chip][b]
            if len(A) != len(B):
                continue
            ch = [abs(y / x - 1) for x, y in zip(A, B) if x > 0 and y > 0]
            if not ch or all(c == 0 for c in ch):
                continue
            print(f'{chip} {a[:16]} -> {b[:16]}  {(P(b) - P(a)).total_seconds() / 3600:6.1f} h apart  median change {100 * S.median(ch):5.1f} %')


if __name__ == '__main__':
    main()
