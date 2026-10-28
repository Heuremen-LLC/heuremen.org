import sys, os, tempfile, shutil, subprocess
sys.path.insert(0, '/Users/heuremen/bones-forge/arxiv_one_readout/flipcheck_run2'); import synth_flip as sf
R = '/Users/heuremen/bones-forge/arxiv_one_readout/'
for sc, seed in (('null', sf.SEED + 100), ('b', sf.SEED + 100 + ord('b'))):
    out = tempfile.mkdtemp(prefix='oldnew_'); sf.make(sc, seed, out)
    for lab, a in (('RUN 1 analyzer (heralded floor)', R + 'flipcheck_run/analyze_flip.py'), ('RUN 2 analyzer (unheralded floor, like for like)', R + 'flipcheck_run2/analyze_flip.py')):
        r = subprocess.run([sys.executable, a, out], capture_output=True, text=True)
        print('#' * 20, 'thermal 1 %% synthetic, scenario %s seed %d, %s' % (sc, seed, lab))
        print('\n'.join(l for l in (r.stdout + r.stderr).splitlines() if 'KNOWN-GOOD' in l or (' F2' in l and ('pooled' in l or 'FINAL' in l)) or 'DECISION' in l or 'Error' in l))
    shutil.rmtree(out)
