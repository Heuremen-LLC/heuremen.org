"""Known-good seat for analyze_flat.py on REAL hardware (CHECK_KNOWNGOOD.md). Runs analyze_flat's own free fit and readings, unchanged,
on the 22 Oct one-neighbour run's data (zz_revival/job_zz_<chip>.json etc.), arms P1/P2 against Z, on its own tau grid."""
import json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.argv = [sys.argv[0], HERE]
import analyze_flat as A
SRC = os.path.join(HERE, '..', 'zz_revival')
class JK(A.Z.Job):
    def __init__(s, chip):
        s.rec = json.load(open(os.path.join(SRC, 'job_zz_%s.json' % chip))); s.z = np.load(os.path.join(SRC, 'raw_zz_%s.npz' % chip))
        s.props = json.load(open(os.path.join(SRC, 'props_zz_%s.json' % chip))); s.trip = [tuple(t) for t in s.rec['info']['triples']]; s.m = len(s.trip)
        s.meta = s.rec['meta']; s.taus = sorted({m['tau_us'] for m in s.meta}); s.edges = {(min(a, c), max(a, c)) for (a, c) in s.rec['edges']}
for chip in ('ibm_fez', 'ibm_phoenix'):
    J = JK(chip); P = A.pairs(J, (0, 1)); ok = np.all(np.isfinite(P), 1)
    print('=' * 90); print('KNOWN-GOOD SEAT %s (22 Oct one-neighbour data, %d of %d pairs fitted)' % (chip, ok.sum(), len(P)))
    for name, f in (('slope of log fitted beat on log calibrated zz', A.slope), ('slope minus 1', lambda x: A.slope(x) - 1),
                    ('rank correlation', A.spear), ('median fitted beat (kHz)', lambda x: A.med(x) / 1e3)):
        v, lo, hi, n = A.boot(P, f); print('  %-48s %+.4g [%+.4g, %+.4g] n=%d -> %s' % (name, v, lo, hi, n, A.oc(lo, hi, n)))
