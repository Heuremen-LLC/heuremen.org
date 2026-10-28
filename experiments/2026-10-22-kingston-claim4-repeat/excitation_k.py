"""Kingston repeat of the claim 4 excitation run. Heuremen LLC, Wayfinder and Bones, 22 Oct 2026. PLAN_EXK.md, sha256 in PLAN_EXK.sha256.
Imports ../deciding_runs/excitation.py (circuit builder, roles, arms) and ../deciding_runs/runs.py (fire path: copies, shuffle seed, transpile,
estimate, submit, save) UNCHANGED, and changes exactly:
  1. shots per circuit: SHOTS below (the original returns 2000; asserted)
  2. chip: ibm_kingston only (the original ran fez, kingston, marrakesh)
  3. output folder: this folder. File names stay job/props/raw_excitation_ibm_kingston.* so ../deciding_runs/analyze_excitation.py reads
     them unchanged (python analyze_excitation.py <this folder>)
  4. the plan hash written into the job file: PLAN_EXK.sha256 (was PLAN_EXCITATION.sha256)
  5. two extra keys in the job file's info record (repeat_of, shots_original); the analysis reads only info['S']
Dry-run by default (builds, transpiles on the live chip, prints ESTIMATED seconds, submits nothing). Usage: python excitation_k.py [--fire]"""
import os, sys, textwrap
HERE = os.path.dirname(os.path.abspath(__file__)); DEC = os.path.normpath(os.path.join(HERE, '..', 'deciding_runs')); sys.path.insert(0, DEC)
import runs, excitation
CHIP = 'ibm_kingston'; SHOTS = 4900; ORIGINAL_JOB = 'dbcnds424mrc73ck0450'
def build_exk(b, chip, shots=None):
    pubs, meta, shots0, info = excitation.build_excitation(b, chip)
    assert shots0 == 2000, 'excitation.py changed: it no longer returns 2000 shots'
    return pubs, meta, SHOTS if shots is None else shots, dict(info, repeat_of=ORIGINAL_JOB, shots_original=shots0)
if __name__ == '__main__':
    runs.BUILD['excitation'] = build_exk
    sys.argv = [sys.argv[0], 'excitation', CHIP] + [a for a in sys.argv[1:] if a == '--fire']
    src = open(runs.__file__).read().split("if __name__ == '__main__':")[1]
    old = "open(HERE + '/PLAN.sha256')"; assert src.count(old) == 1; src = src.replace(old, "open(HERE + '/PLAN_EXK.sha256')")
    g = dict(runs.__dict__); g['HERE'] = HERE   # every output path in the fire path is HERE + '/...': this folder, never ../deciding_runs
    exec(textwrap.dedent(src), g)
