"""Post-hoc reanalysis (after the data were seen): random-excitation prediction for the readout coherence excess.
A measured qubit left in 1 with probability p applies its phase only in those shots, costing about p*(1 - cos(phi)) of coherence,
where phi is the phase measured in the idle-in-1 arm over 8 steps (shared equally among a spectator's coupled measured qubits).
p is measured from the recorded mid-circuit readouts in the readout-in-0 arm."""
import json, os, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); rng = np.random.default_rng(9)
R = os.path.join(HERE, '..', '2026-10-04-three-positions', 'data', 'reads', 'countwide_roles_%s.json')
def bits(arr, k):
    nby = arr.shape[1]; return ((arr[:, nby - 1 - k // 8] >> (k % 8)) & 1).astype(int)
def p1(arr, n): return np.array([bits(arr, k).mean() for k in range(n)])
def med(x): b = [np.median(rng.choice(x, len(x))) for _ in range(2000)]; return np.median(x), np.percentile(b, [2.5, 97.5])
for chip in ['ibm_kingston', 'ibm_marrakesh']:
    rec = json.load(open(os.path.join(HERE, 'job_%s.json' % chip))); z = np.load(os.path.join(HERE, 'raw_%s.npz' % chip)); E = rec['E']; S = rec['S']; meta = rec['meta']
    nb = {}
    for a, b in json.load(open(R % chip))['edges']: nb.setdefault(a, set()).add(b); nb.setdefault(b, set()).add(a)
    P = {}; pE = []
    for i, m in enumerate(meta):
        P[(m['mode'], m['N'], m['quad'])] = p1(z['fin_%d' % i], len(S))
        if m['mode'] == 'M0':
            for k in range(m['N']):
                key = 'm%d_%d' % (k, i)
                if key in z.files: pE.append(p1(z[key], len(E)))
    pE = np.mean(pE, axis=0); Ei = {q: j for j, q in enumerate(E)}
    V = lambda mo, N: (1 - 2 * P[(mo, N, 'X')]) + 1j * (1 - 2 * P[(mo, N, 'Y')])
    coh = lambda mo, N: np.abs(V(mo, N))
    c0 = coh('D0', 0); Xex = (c0 - coh('M0', 8)) - (c0 - coh('D0', 8))
    phi = np.angle(V('D1', 8) / V('D0', 8))
    coupled = np.array([len(nb.get(q, set()) & set(E)) > 0 for q in S])
    nE = np.array([len(nb.get(q, set()) & set(E)) for q in S])
    pnb = np.array([np.mean([pE[Ei[e]] for e in nb.get(q, set()) & set(E)]) if c else 0.0 for q, c in zip(S, coupled)])
    pred = np.array([nE[j] * pnb[j] * (1 - np.cos(phi[j] / nE[j])) if coupled[j] else 0.0 for j in range(len(S))]) * coh('D0', 8)
    print('%s: measured-qubit excitation rate %.4f' % (chip, pE.mean()))
    for lab, mask in (('next to a measured qubit', coupled), ('no measured neighbour', ~coupled)):
        mx, cx = med(Xex[mask]); mp, _ = med(pred[mask]); md, cd = med((Xex - pred)[mask])
        print('   %-25s n=%3d  excess %+.4f [%+.4f,%+.4f]  random-excitation prediction %+.4f  difference %+.4f [%+.4f,%+.4f]' % (lab, mask.sum(), mx, *cx, mp, md, *cd))
