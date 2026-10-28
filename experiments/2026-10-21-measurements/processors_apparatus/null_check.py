"""Is the share of qubit-copies below their own single-exponential extrapolation more than noise would give?
For each qubit, the three copies (a, b, c) in the same job give the noise at 30 us and 300 us. A pure exponential
through the mean 300 us level is drawn 20 times per qubit with that noise, and the same 'below' test is applied.
Standard library only. Reads data/time_runs_ruler_copies.json."""
import json, random, statistics as st
d = json.load(open('data/time_runs_ruler_copies.json')); rng = random.Random(1)
lo = [1, 0]
for run in d['runs']:
    obs, nul = [], []
    for r in run['rows']:
        cps = [r[1:5], r[5:9], r[9:13]]
        obs += [c[0] < c[2] ** 0.1 for c in cps if c[2] > 0]
        s30 = [c[0] for c in cps]; s300 = [c[2] for c in cps]
        if min(s300) <= 0: continue
        m300 = st.mean(s300); t30 = m300 ** 0.1
        for _ in range(20):
            b = m300 + rng.gauss(0, st.stdev(s300))
            if b > 0: nul.append(t30 + rng.gauss(0, st.stdev(s30)) < b ** 0.1)
    o, n = sum(obs) / len(obs), sum(nul) / len(nul)
    print('%-15s run %4d  below own extrapolation: observed %.1f%%   pure-exponential null %.1f%%' % (run['chip'], run['run'], 100 * o, 100 * n))
