"""The 4 October spaced X-gate sweep on ibm_pittsburgh with all three jobs (runs 361, 362, 363; same circuits) against the first job alone,
with verify_xtrain.py's own functions. Defined after the data; the paper first used run 361 only."""
import os
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'verify_xtrain.py')).read().split("def main():")[0]
exec(src)
for runs in ((361,), (361, 362, 363)):
    rs = sweep_chip('ibm_pittsburgh', runs)
    print('pittsburgh runs', runs)
    for t in (5, 10, 20, 30, 150):
        m1, w1 = sw(rs, 'OFF', t, 'ONE'); m2, w2 = sw(rs, 'ON', t, 'ONE'); c1, _ = sw(rs, 'OFF', t, 'CURVE')
        print('  %3d us  start0 vs one exponential %+.4f +-%.4f | start1 vs one exponential %+.4f | start0 vs own curve %+.4f' % (t, m1, w1, m2, c1))
