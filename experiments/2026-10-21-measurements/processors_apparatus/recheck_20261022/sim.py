import json, numpy as np, sys, itertools
D='/Users/heuremen/bones-forge/heuremen.org/experiments/2026-10-21-measurements/processors_apparatus/data/time_runs_ruler_copies.json'
runs=json.load(open(D))['runs']
T=np.array([30,100,300,600.]); EXP=30/300
rng=np.random.default_rng(7)

def fitT1(Lm):
    g=np.linspace(40,900,431); best=None
    err=((np.exp(-T[None,:]/g[:,None])-Lm[None,:])**2).sum(1)
    return g[err.argmin()]

def chain(pop, n, sep, fd, epp, e0, R, shared=True, mode='iq'):
    """pop: (Q,4) true excited-state population (relative to ON ref), returns levels (R,Q,3,4) through IQ chain.
    1-D projection; OFF centre 0, ON centre sep, width 1; fd = fraction of excited shots that relax mid-readout (land between).
    ON ref = prepared 1 with prep error epp; OFF ref = prepared 0 with residual e0. level=(mean-offmean)/(onmean-offmean)."""
    Q=pop.shape[0]
    def meanproj(p, shape):
        # p: excited probability array broadcast to shape; return mean projection of n shots (normal approx of the sum)
        k=rng.binomial(n,np.broadcast_to(p,shape))            # excited shots
        kd=rng.binomial(k,fd); kn=k-kd
        s=(kn*sep + kd*sep*0.5 + rng.normal(0,np.sqrt(np.maximum(kd,0)*(sep**2/12)))) 
        s+= rng.normal(0,np.sqrt(n),shape)                    # unit-width Gaussian noise on every shot
        return s/n
    if mode=='bits':
        # classified bits, threshold midway; raw fraction of 1s (no normalisation)
        def frac(p,shape):
            k=rng.binomial(n,np.broadcast_to(p,shape)); kd=rng.binomial(k,fd); kn=k-kd
            # P(read 1| excited, no decay)=Phi(sep/2); decayed: uniform centre -> ~0.5 ; ground: Phi(-sep/2)
            from math import erf,sqrt
            ph=0.5*(1+erf((sep/2)/sqrt(2)))
            r1=rng.binomial(kn,ph)+rng.binomial(kd,0.5)+rng.binomial(n-k,1-ph)
            return r1/n
        meas=frac
    else: meas=meanproj
    shp=(R,Q,3,4)
    pt=pop*(1-epp)+0   # idle state starts from ON-prepared state (1-epp)
    # idle state population at t: decays from (1-epp) toward e0 thermal
    pt=e0+(1-epp-e0)*pop
    x=meas(pt[None,:,None,:],shp)
    if shared:
        on=meas(np.full((Q,),1-epp),(R,Q))[:,:,None,None]; off=meas(np.full((Q,),e0),(R,Q))[:,:,None,None]
    else:
        on=meas(np.full((Q,),1-epp)[None,:,None,None],(R,Q,3,1)); off=meas(np.full((Q,),e0)[None,:,None,None],(R,Q,3,1))
    if mode=='bits': return x
    return (x-off)/(on-off)

def below(L):  # L (...,3,4): copy-wise test used in the claim
    c30=L[...,0]; c300=L[...,2]
    ok=c300>0
    return (c30<np.where(ok,c300,1)**EXP), ok

def calib_n(Lobs, sep, fd, epp, e0, pop):
    # match mean within-qubit copy variance pooled over the 4 times
    target=(Lobs.std(1,ddof=1)**2).mean(0)  # (4,)
    best=None
    for n in range(100,1500,100):
        L=chain(pop,n,sep,fd,epp,e0,8,shared=True)
        v=(L.std(2,ddof=1)**2).mean((0,1))
        e=((np.log(v)-np.log(target))**2).sum()
        if best is None or e<best[0]: best=(e,n)
    return best[1]

def published_null(Lc, r):
    # exact logic of null_check.py on a (Q,3,4) array
    nul=[]
    for q in range(Lc.shape[0]):
        s30=Lc[q,:,0]; s300=Lc[q,:,2]
        if s300.min()<=0: continue
        m=s300.mean(); t=m**EXP
        for _ in range(20):
            b=m+r.gauss(0,np.std(s300,ddof=1))
            if b>0: nul.append(t+r.gauss(0,np.std(s30,ddof=1))<b**EXP)
    return np.mean(nul)
