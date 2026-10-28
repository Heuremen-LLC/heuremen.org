"""Analysis exactly as PREREG.md (sha256 in PREREG.sha256)."""
import json, glob, os, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); rng = np.random.default_rng(7)
ROLES = os.path.join(HERE, '..', '2026-10-04-three-positions', 'data', 'reads', 'countwide_roles_%s.json')
def bits(arr, k):
    nby = arr.shape[1]; return ((arr[:, nby - 1 - k // 8] >> (k % 8)) & 1)
def p1(arr, nb): return np.array([bits(arr, k).mean() for k in range(nb)])
def boot(x, n=2000):
    x = np.asarray(x); b = [np.median(rng.choice(x, len(x))) for _ in range(n)]; return np.median(x), np.percentile(b, [2.5, 97.5])
def state(lo, hi): return 'ABOVE' if lo > 0 else 'BELOW' if hi < 0 else 'CONSISTENT'
for jf in sorted(glob.glob(HERE + '/job_*.json')):
    rec = json.load(open(jf)); chip = rec['chip']; raw = HERE + '/raw_%s.npz' % chip
    if not os.path.exists(raw): print(chip, 'NOTHING_RETURNED'); continue
    z = np.load(raw); E = rec['E']; S = rec['S']; meta = rec['meta']
    roles = json.load(open(ROLES % chip)); nb = {}
    for a, b in roles['edges']: nb.setdefault(a, set()).add(b); nb.setdefault(b, set()).add(a)
    P = {}; pE = []
    for i, m in enumerate(meta):
        P[(m['mode'], m['N'], m['quad'])] = p1(z['fin_%d' % i], len(S))
        if m['mode'] == 'M0':
            for k in range(m['N']):
                key = 'm%d_%d' % (k, i)
                if key in z.files: pE.append(p1(z[key], len(E)))
    pE = np.mean(pE, axis=0)  # per-emitter excitation probability in M0, averaged over steps and pubs
    coh = lambda mo, N: np.abs((1 - 2 * P[(mo, N, 'X')]) + 1j * (1 - 2 * P[(mo, N, 'Y')]))
    c0 = coh('D0', 0)
    L = {mo: c0 - coh(mo, 8) for mo in ('D0', 'M0', 'D1', 'M1')}
    X = L['M0'] - L['D0']; Y = L['D1'] - L['D0']
    Eidx = {q: j for j, q in enumerate(E)}
    coupled = np.array([len(nb.get(q, set()) & set(E)) > 0 for q in S])
    pnb = np.array([np.mean([pE[Eidx[e]] for e in nb.get(q, set()) & set(E)]) if c else pE.mean() for q, c in zip(S, coupled)])
    Xp = pnb * Y; D = X - Xp
    print('%s: E excitation in M0 (mean over emitters, steps) p = %.4f  (max emitter %.3f); spectators %d (%d coupled to an E)' % (chip, pE.mean(), pE.max(), len(S), coupled.sum()))
    for name, mask in (('all', np.ones(len(S), bool)), ('coupled to an E', coupled), ('not coupled', ~coupled)):
        if mask.sum() < 3: print('   [%s] n=%d CANNOT_TELL' % (name, mask.sum())); continue
        mx, cx = boot(X[mask]); my, cy = boot(Y[mask]); md, cd = boot(D[mask]); mp, cp = boot(Xp[mask])
        verdict = 'CANNOT_TELL (no readout excess to explain)' if cx[0] <= 0 <= cx[1] else state(*cd)
        print('   [%-15s] n=%3d  X=%+.4f [%+.4f,%+.4f]  Y=%+.4f [%+.4f,%+.4f]  X_pred=%+.4f  X-X_pred=%+.4f [%+.4f,%+.4f] -> %s' % (name, mask.sum(), mx, *cx, my, *cy, mp, md, *cd, verdict))
    print('   (also printed) M1 readout excess over D1: %+.4f' % np.median(L['M1'] - L['D1']))
