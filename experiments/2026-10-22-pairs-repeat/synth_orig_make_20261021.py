import json, numpy as np, random
rng = np.random.default_rng(11); random.seed(11)
n = 156; E = list(range(39)); S = list(range(39, 156)); P = [[39 + 2 * j, 40 + 2 * j] for j in range(48)]; SG = list(range(135, 156)); K1 = [0]
edges = [[e, s] for e in E for s in random.sample(S[:96], 2)]
json.dump({'gates': [{'gate': 'cz', 'qubits': ed} for ed in edges]}, open('pairs_backend_properties_synthetic.json', 'w'))
adj = {}
for a, b in edges: adj.setdefault(a, set()).add(b); adj.setdefault(b, set()).add(a)
nbP = [[E.index(e) for e in (adj.get(a, set()) | adj.get(b, set())) if e in E] for a, b in P]
mu0, mu1, mu2, sg = 1 + 0j, -1 + 0j, 1 + 1j, 0.2; SH = 800
iq = lambda centre, size: centre + sg * (rng.standard_normal(size) + 1j * rng.standard_normal(size))
for copy in (0, 1):
    meta, data = [], {}
    def add(m, regs):
        i = len(meta); meta.append(dict(m, copy=copy)); data.update({'%s_%d' % (k, i): v.astype(np.complex64) for k, v in regs.items()})
    for st in '01':
        bits = np.full((SH, n), int(st)); add(dict(arm='CAL', state=st), {'c': np.where(bits == 1, iq(mu1, (SH, n)), iq(mu0, (SH, n)))})
    for N, mode, K in [(0, 'D', 39), (2, 'D', 39), (2, 'M', 39), (4, 'D', 39), (4, 'M', 39), (4, 'M', 1)]:
        for st in ('ZZ', 'XX', 'YY', 'XY', 'YX'):
            Kn = 39 if K == 39 else 1
            leak = (rng.random((SH, Kn)) < 0.03 * N / 4) if (mode == 'M' and K == 39) else np.zeros((SH, Kn), bool)
            base_c = max(0.0, 0.8 - 0.1 * N); loss = 0.03 if (mode == 'M' and K == 39 and N > 0) else 0.0
            b = np.zeros((SH, 2 * 48 + 21), int)
            for j in range(48):
                c = np.full(SH, max(0.0, base_c - loss))
                if K == 39 and nbP[j]: c[leak[:, nbP[j]].any(1)] = 0.0
                b1 = rng.integers(0, 2, SH)
                same = rng.random(SH) < (0.95 if st == 'ZZ' else (1 + {'XX': c, 'YY': -c, 'XY': 0 * c, 'YX': 0 * c}[st]) / 2)
                b[:, 2 * j] = b1; b[:, 2 * j + 1] = np.where(same, b1, 1 - b1)
            b[:, 96:] = rng.integers(0, 2, (SH, 21))
            regs = {'fin': np.where(b == 1, iq(mu1, b.shape), iq(mu0, b.shape))}
            for k in range(N if mode == 'M' else 0): regs['m%d' % k] = np.where(leak, iq(mu2, (SH, Kn)), iq(mu0, (SH, Kn)))
            add(dict(arm='P', mode=mode, N=N, K=K, setting=st), regs)
    json.dump(dict(chip='synthetic', copy=copy, job_id='none', meta=meta, E=E, S=S, pairs=P, singles=SG, K1=K1), open('pairsjob_synthetic_copy%d.json' % copy, 'w'))
    np.savez_compressed('pairsraw_synthetic_copy%d.npz' % copy, **data)
print('ok')
