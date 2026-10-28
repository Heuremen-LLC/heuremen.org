"""Synthetic test of control_run/analyze_pairs_iq.py (UNCHANGED) in the REAL pairs schema, for the repeat (pairs_repeat/PREREG_PAIRS.md).
Heuremen LLC, Wayfinder and Bones, 22 October 2026.
The 21 Oct test (synth_orig_20261021/make.py) used a toy chip (800 shots, loss 0.03) and had no null; this one sets the floor at the real scale.
Per chip (ibm_fez, ibm_kingston, ibm_marrakesh): roles E, S, pairs, singles, K1 and the meta schema are the REAL 21 Oct ones
(control_run/pairsjob_<chip>_copy<c>.json, the 21 Oct pub order per copy); cz wiring is the REAL 21 Oct properties file; every qubit's 0/1 IQ
centre and spread come from the REAL 21 Oct calibration pubs (pairsraw_<chip>_copy0.npz, c_ register). Registers, shapes and dtype as the
runner writes them: c_i (2000 x 156) for CAL, fin_i (4000 x 117: pair members 2j, 2j+1, then singles), m<k>_i (4000 x 39 for K=39, 4000 x 1
for K=1; all zeros in D mode, as on hardware). Files are written with np.savez (uncompressed; np.load reads both the same way).
Common to all scenarios: per pair idle XX/YY correlation c0 ~ U(0.85, 0.97) x exp(-0.004 N), ZZ agreement 0.97; singles random; measured
qubits read in their 0 cluster on mid readouts; in each M-mode mid readout each measured qubit lands VISIBLY in its 1 cluster with
probability 0.25 % (an excitation the clean-shot selection can see).
Scenarios (chip name in the files is synth-<chip>, so nothing can be mistaken for hardware):
  a     PLANTED: wired pairs (>= 1 measured cz neighbour) lose 0.009 of c at N=4 K=39 and 0.0045 at N=2 K=39, 0.001 at N=4 K=1, in EVERY shot;
        in shots where a wired neighbour was visibly excited the pair loses a further 0.5 (so all-shot R3 > clean-shot R1; R1 should read ~0.009)
  leak  loss ONLY in shots with a visibly excited wired neighbour (0.5 there, 0 elsewhere): R1 should NOT be ABOVE, R3 should be
  null  nothing planted (no loss, visible excitations carry no loss)
Usage: synth_pairs.py SCENARIO SEED [chip ...]   -> writes synth_<scenario>_<seed>.txt here; raw files go to a temp folder and are deleted."""
import os, sys, json, subprocess, tempfile, shutil
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); CTRL = os.path.join(os.path.dirname(HERE), 'control_run'); ANALYSIS = os.path.join(CTRL, 'analyze_pairs_iq.py')
PY = sys.executable; SC = sys.argv[1]; SEED = int(sys.argv[2]); CHIPS = sys.argv[3:] or ['ibm_fez', 'ibm_kingston', 'ibm_marrakesh']
assert SC in ('a', 'leak', 'null')
PEXC = 0.0025; KICK = 0.5; LOSS = {(4, 39): 0.009, (2, 39): 0.0045, (4, 1): 0.001}

def make(chip, out, rng):
    props = json.load(open(os.path.join(CTRL, 'pairs_backend_properties_%s.json' % chip)))
    json.dump({'gates': [g for g in props['gates'] if g['gate'] == 'cz']}, open(os.path.join(out, 'pairs_backend_properties_synth-%s.json' % chip), 'w'))
    adj = {}
    for g in props['gates']:
        if g['gate'] == 'cz': a, b = g['qubits']; adj.setdefault(a, set()).add(b); adj.setdefault(b, set()).add(a)
    zr = np.load(os.path.join(CTRL, 'pairsraw_%s_copy0.npz' % chip)); r0 = json.load(open(os.path.join(CTRL, 'pairsjob_%s_copy0.json' % chip)))
    C0 = np.concatenate([zr['c_%d' % i] for i, m in enumerate(r0['meta']) if m['arm'] == 'CAL' and m['state'] == '0'])
    C1 = np.concatenate([zr['c_%d' % i] for i, m in enumerate(r0['meta']) if m['arm'] == 'CAL' and m['state'] == '1'])
    MU0 = np.median(C0.real, 0) + 1j * np.median(C0.imag, 0); MU1 = np.median(C1.real, 0) + 1j * np.median(C1.imag, 0)
    SIG = np.median(np.abs(C0 - MU0), 0) / np.sqrt(2 * np.log(2))   # per-quadrature sigma from the median radius of a 2-D Gaussian
    E, S, P, SG, K1 = r0['E'], r0['S'], r0['pairs'], r0['singles'], r0['K1']; n = len(E) + len(S); nP = len(P)
    Epos = {e: i for i, e in enumerate(E)}; nb = [sorted(Epos[e] for e in (adj.get(a, set()) | adj.get(b, set())) if e in Epos) for a, b in P]
    wired = np.array([len(x) > 0 for x in nb]); c0 = rng.uniform(0.85, 0.97, nP); pz = 0.97
    def iq(q, bits):   # q: qubit numbers per column; bits: shots x cols in {0,1}
        q = np.asarray(q); mu = np.where(bits == 1, MU1[q], MU0[q])
        return (mu + SIG[q] * (rng.standard_normal(bits.shape) + 1j * rng.standard_normal(bits.shape))).astype(np.complex64)
    for copy in (0, 1):
        rec = json.load(open(os.path.join(CTRL, 'pairsjob_%s_copy%d.json' % (chip, copy)))); data = {}
        for i, m in enumerate(rec['meta']):
            if m['arm'] == 'CAL':
                data['c_%d' % i] = iq(np.arange(n), np.full((2000, n), int(m['state']))); continue
            mode, N, K, st = m['mode'], m['N'], m['K'], m['setting']; sh = 4000
            Kset = E if K == 39 else K1; exc = np.zeros((sh, len(Kset)), bool)
            for k in range(N):
                if mode == 'M':
                    ek = rng.random((sh, len(Kset))) < PEXC; exc |= ek; data['m%d_%d' % (k, i)] = iq(Kset, ek.astype(int))
                else: data['m%d_%d' % (k, i)] = np.zeros((sh, len(Kset)), np.complex64)
            kpos = {Epos[e]: j for j, e in enumerate(Kset)}
            bits = np.zeros((sh, 2 * nP + len(SG)), int)
            for j in range(nP):
                c = np.full(sh, c0[j] * np.exp(-0.004 * N))
                if mode == 'M' and wired[j]:
                    hit = [kpos[x] for x in nb[j] if x in kpos]; bad = exc[:, hit].any(1) if hit else np.zeros(sh, bool)
                    if SC == 'a': c -= LOSS.get((N, K), 0.0); c[bad] -= KICK
                    elif SC == 'leak': c[bad] -= KICK
                c = np.clip(c, -1, 1)
                b1 = rng.integers(0, 2, sh)
                psame = {'ZZ': np.full(sh, pz), 'XX': (1 + c) / 2, 'YY': (1 - c) / 2, 'XY': np.full(sh, 0.5), 'YX': np.full(sh, 0.5)}[st]
                same = rng.random(sh) < psame; bits[:, 2 * j] = b1; bits[:, 2 * j + 1] = np.where(same, b1, 1 - b1)
            bits[:, 2 * nP:] = rng.integers(0, 2, (sh, len(SG)))
            data['fin_%d' % i] = iq([q for pr in P for q in pr] + list(SG), bits)
        rec = dict(rec, chip='synth-%s' % chip, job_id='SYNTH-%s-%d' % (SC, SEED), prereg_sha256='synthetic')
        json.dump(rec, open(os.path.join(out, 'pairsjob_synth-%s_copy%d.json' % (chip, copy)), 'w'))
        np.savez(os.path.join(out, 'pairsraw_synth-%s_copy%d.npz' % (chip, copy)), **data)
    return int(wired.sum())

if __name__ == '__main__':
    rng = np.random.default_rng(SEED); tmp = tempfile.mkdtemp(prefix='synth_pairs_')
    try:
        nw = {chip: make(chip, tmp, rng) for chip in CHIPS}
        res = subprocess.run([PY, ANALYSIS, tmp], capture_output=True, text=True)
        hdr = ('SYNTHETIC scenario %s seed %d (synth_pairs.py; analysis = control_run/analyze_pairs_iq.py unchanged). Wired pairs per chip: %s\n'
               'Planted (scenario a): clean-shot c loss N=4 0.009, N=2 0.0045 (R1 reads it times the readout contrast, ~0.9); K=1 0.001; R2 and R3 also carry the\n'
               'visible-excitation kicks (all shots). leak: loss only in visibly-excited-neighbour shots. null: nothing planted.\n'
               % (SC, SEED, nw))
        open(os.path.join(HERE, 'synth_%s_%d.txt' % (SC, SEED)), 'w').write(hdr + res.stdout + res.stderr)
        print(hdr + res.stdout + res.stderr)
    finally: shutil.rmtree(tmp)
