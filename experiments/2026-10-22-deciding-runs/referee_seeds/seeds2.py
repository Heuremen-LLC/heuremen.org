import sys,io,contextlib,re,pickle,time,collections
D='/Users/heuremen/bones-forge/arxiv_one_readout/deciding_runs/'
sys.path.insert(0,D); F=D+'only_deciding'; sys.argv=['x',F]
import numpy as np, analyze as A
chips=['ibm_fez','ibm_kingston','ibm_marrakesh','ibm_phoenix']
runs=['order','neighbour','gates']
seeds=int(sys.argv_n) if False else 40
tab=collections.defaultdict(collections.Counter)
comb=collections.defaultdict(list)
for s in range(1,seeds+1):
    A.RNG=np.random.default_rng(s); A.LOG.clear(); buf=io.StringIO()
    with contextlib.redirect_stdout(buf):
        for run in runs:
            for chip in chips:
                try: J=A.Job(run,chip,F)
                except FileNotFoundError: continue
                J.props_n=J.rec.get('n_qubits',156); A.CUR[0]=chip; A.RUNS[run](J)
    for lab,chip,p in A.LOG:
        tab[(run_of:=lab,chip)][ 'CANNOT TELL(shift)' if p is None else p[4]]+=1
    by=collections.defaultdict(list)
    for lab,chip,p in A.LOG: by[lab].append(p)
    # combined outcome
    for lab,rows in by.items():
        ok=[p for p in rows if p is not None and np.isfinite(p[1])]
        if len(ok)<2: continue
        v=np.array([p[0] for p in ok]); sg=np.array([(p[2]-p[1])/3.92 for p in ok]); w=1/sg**2; m=(w*v).sum()/w.sum(); se=1/np.sqrt(w.sum())
        tab[(lab,'COMBINED')][A.outcome(m-1.96*se,m+1.96*se,99)]+=1
pickle.dump(dict(tab),open('tab.pkl','wb'))
for k,c in sorted(tab.items()):
    if len(c)>1: print(k,dict(c))
print('done')
