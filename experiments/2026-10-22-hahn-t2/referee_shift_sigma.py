"""Outside referee, 23 Oct 2026 (defined after the data): run-1 to run-2 shift of the planned-fit median under one unchanged record, against shot noise alone (binomial redraw of both runs, 300 draws). Kingston and marrakesh."""
import json, numpy as np
D='/Users/heuremen/bones-forge/heuremen.org/experiments/2026-10-22-hahn-t2/'
T=np.array([0,10,20,40,60,90,120,160,220,300,400.]); G=np.geomspace(1,2000,600); E=np.exp(-T[None,:]/G[:,None])
def p1(tag):
    job=json.load(open(D+f'job_hahn_ibm_{tag}.json')); z=np.load(D+f'raw_hahn_ibm_{tag}.npz'); acc={}
    for i,m in enumerate(job['meta']):
        S=job['sets'][m['set']]; b=z[f'c_{i}'][:,:len(S)]
        for j,q in enumerate(S): acc.setdefault(q,{}).setdefault(m['tau_us'],[]).append(b[:,j].mean())
    return {q:np.array([np.mean(acc[q][t]) for t in T]) for q in acc}
def fitmany(Y):  # Y: (n,11) fixed fit
    yy=Y-.5; A=(yy@E.T)/(E*E).sum(1); r=((yy[:,None,:]-A[:,:,None]*E[None])**2).sum(2); k=r.argmin(1); return G[k],A[np.arange(len(k)),k]
def med(P,qs,st):
    Y=np.array([P[q] for q in qs]); T2,A=fitmany(Y); m=(A>0.2); return np.median(T2[m]/st[m])
rng=np.random.default_rng(3)
for chip,a,b,pr in [('kingston','kingston_run1','kingston','kingston_run1'),('marrakesh','marrakesh_run1','marrakesh','marrakesh_run1')]:
    props=json.load(open(D+f'props_hahn_ibm_{pr}.json'))
    stated={q:[x['value'] for x in ps if x['name']=='T2'][0] for q,ps in enumerate(props['qubits']) if any(x['name']=='T2' for x in ps)}
    A_,B_=p1(a),p1(b); qs=sorted(set(A_)&set(B_)&set(stated)); st=np.array([stated[q] for q in qs])
    obs=med(B_,qs,st)-med(A_,qs,st)
    N=4000  # 2 copies x 2000 shots per tau
    d=[]
    for _ in range(300):
        Ar={q:rng.binomial(N,np.clip(A_[q],0,1))/N for q in qs}; Br={q:rng.binomial(N,np.clip(B_[q],0,1))/N for q in qs}
        d.append(med(Br,qs,st)-med(Ar,qs,st))
    d=np.array(d); print(chip,'obs shift %.3f'%obs,'shot-noise sd of shift %.4f'%d.std(),'z %.1f'%(obs/d.std()))
