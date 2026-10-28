"""EXPLORATORY -- NOT in PREREG.md (written after the data). Splits spectators into neighbours of a measured qubit vs far qubits."""
import json, os, numpy as np
rng = np.random.default_rng(11); HERE = os.path.dirname(os.path.abspath(__file__))
ROLES = os.path.join(HERE, 'roles_%s.json')
def p1(arr, nb):
    nby = arr.shape[1]; return np.array([((arr[:, nby - 1 - k // 8] >> (k % 8)) & 1).mean() for k in range(nb)])
def slope(ns, ph): return np.polyfit(ns, np.unwrap(ph), 1)[0]
for chip in ('ibm_fez', 'ibm_kingston', 'ibm_pittsburgh', 'ibm_marrakesh'):
    jf = os.path.join(HERE, 'job_%s.json' % chip); rf = os.path.join(HERE, 'raw_%s.npz' % chip)
    if not (os.path.exists(jf) and os.path.exists(rf)): continue
    rec = json.load(open(jf)); z = np.load(rf); S = rec['S']; meta = rec['meta']
    nb = set(int(k) for k in json.load(open(ROLES % chip))['witness_of'].keys())
    P = {(m['mode'], m['N'], m['quad']): p1(z['fin_%d' % i], len(S)) for i, m in enumerate(meta)}
    NS = sorted({m['N'] for m in meta}); s = {}
    for mode in ('D', 'M0', 'M1'):
        ph = np.array([np.arctan2(2 * P[(mode, N, 'Y')] - 1, 2 * P[(mode, N, 'X')] - 1) for N in NS])
        s[mode] = np.array([slope(np.array(NS, float), ph[:, k]) for k in range(len(S))])
    keep = np.abs(np.degrees(s['D'])) >= 1.0
    for lab, sel in (('neighbour', np.array([q in nb for q in S])), ('far', np.array([q not in nb for q in S]))):
        k = keep & sel; R0 = s['M0'][k] / s['D'][k]; d = (s['M1'][k] - s['M0'][k]) / s['D'][k]
        b = [np.median(rng.choice(R0, len(R0))) for _ in range(2000)]; bd = [np.median(rng.choice(d, len(d))) for _ in range(2000)]
        print('%-14s %-9s n=%3d  R0 %.3f 95%% [%.3f,%.3f]   R1-R0 %+.3f 95%% [%+.3f,%+.3f]' % (chip, lab, k.sum(), np.median(R0), *np.percentile(b, [2.5, 97.5]), np.median(d), *np.percentile(bd, [2.5, 97.5])))
