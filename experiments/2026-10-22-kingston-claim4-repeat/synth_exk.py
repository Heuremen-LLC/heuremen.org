"""Synthetic test for the kingston repeat. Heuremen LLC, Wayfinder and Bones, 22 Oct 2026. Runs ../deciding_runs/synth_excitation.py's own
source with exactly these substitutions: FakeFez -> FakeKingston and roles ibm_fez -> ibm_kingston (117 spectators, as on the chip); the builder
call goes through excitation_k.build_exk at SHOTS; the two planted cases become CASES (G = 0 and G = +0.5 deg/step); the noise seed 9 becomes
SEED (1..6); outputs go to a temporary folder. The analysis is analyze_excitation.excitation, unchanged. Usage: python synth_exk.py [shots]
NOTE: the synthetic gives every spectator the same G, so its interval is the SHOT-NOISE-ONLY width (a floor). Real spectators differ."""
import os, sys, tempfile, shutil
HERE = os.path.dirname(os.path.abspath(__file__)); DEC = os.path.normpath(os.path.join(HERE, '..', 'deciding_runs')); sys.path.insert(0, HERE); sys.path.insert(0, DEC)
import numpy as np
import excitation_k as EXK, analyze as A
SHOTS = int(sys.argv[1]) if len(sys.argv) > 1 else EXK.SHOTS; SEEDS = range(1, 7); CASES = (('null', 0.0), ('planted', 0.5))
src = open(os.path.join(DEC, 'synth_excitation.py')).read()
subs = [("for label, G in (('planted', -0.6), ('null', 0.0)):", 'for label, G in CASES:'),
        ('rng = np.random.default_rng(9)', 'rng = np.random.default_rng(SEED)'),
        ("OUT = os.path.join(HERE, 'synth_exc_' + label)", "OUT = os.path.join(OUTROOT, label)"),
        ("excitation.build_excitation(b, 'ibm_fez')", "EXK.build_exk(b, 'ibm_kingston', SHOTS)"),
        ('FakeFez', 'FakeKingston'), ("'ibm_fez'", "'ibm_kingston'"), ('raw_excitation_ibm_fez', 'raw_excitation_ibm_kingston'),
        ('props_excitation_ibm_fez', 'props_excitation_ibm_kingston'), ('job_excitation_ibm_fez', 'job_excitation_ibm_kingston')]
for a, b_ in subs:
    assert a in src, a; src = src.replace(a, b_)
src = src.replace("print('#', label, 'G planted', G);", "A.CUR[0] = '%s seed %d' % (label, SEED); print('#', label, 'G planted', G, 'seed', SEED, 'shots', shots);")
for SEED in SEEDS:
    OUTROOT = tempfile.mkdtemp()
    g = dict(__file__=os.path.join(DEC, 'synth_excitation.py'), __name__='synth_exk', CASES=CASES, SEED=SEED, OUTROOT=OUTROOT, SHOTS=SHOTS, EXK=EXK)
    exec(compile(src, 'synth_excitation.py (substituted)', 'exec'), g); shutil.rmtree(OUTROOT)
print('=' * 100); print('SUMMARY, %d shots, seeds %s' % (SHOTS, list(SEEDS)))
for case, G in CASES:
    for lab in sorted({l for l, _, _ in A.LOG}):
        rows = [(c, p) for l, c, p in A.LOG if l == lab and c.startswith(case + ' ')]
        ok = [p for _, p in rows if p is not None]; flags = len(rows) - len(ok)
        if not ok: print('  %-8s %-50s copies shifted in all %d seeds' % (case, lab, len(rows))); continue
        v = np.array([p[0] for p in ok]); w = np.array([p[2] - p[1] for p in ok]); out = [p[4] for p in ok]
        print('  %-8s %-50s value mean %+.4f (min %+.4f max %+.4f) | width mean %.4f (min %.4f max %.4f) | outcomes %s | copy flags %d/%d'
              % (case, lab[:50], v.mean(), v.min(), v.max(), w.mean(), w.min(), w.max(), {o: out.count(o) for o in sorted(set(out))}, flags, len(rows)))
