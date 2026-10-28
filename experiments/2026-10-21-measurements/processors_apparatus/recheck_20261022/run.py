from sim import *
import random, statistics as st
V={ 'main':dict(sep=6,fd=.03,epp=.01,e0=.005,shared=True,floor=True),
 'sep4':dict(sep=4,fd=.03,epp=.01,e0=.005,shared=True,floor=True),
 'sep9':dict(sep=9,fd=.03,epp=.01,e0=.005,shared=True,floor=True),
 'nodecay':dict(sep=6,fd=0,epp=.01,e0=.005,shared=True,floor=True),
 'decay10':dict(sep=6,fd=.10,epp=.01,e0=.005,shared=True,floor=True),
 'prep5':dict(sep=6,fd=.03,epp=.05,e0=.005,shared=True,floor=True),
 'unshared':dict(sep=6,fd=.03,epp=.01,e0=.005,shared=False,floor=True),
 'nofloor':dict(sep=6,fd=.03,epp=.01,e0=.005,shared=True,floor=False),
 'worst':dict(sep=4,fd=.10,epp=.05,e0=.01,shared=True,floor=True)}
R=300; B=2000
out=[]; pr=random.Random(1)
for run in runs:
    A=np.array(run['rows']); Q=len(A)
    Lobs=np.stack([A[:,1:5],A[:,5:9],A[:,9:13]],1)    # (Q,3,4)
    idle0=A[:,13].mean()
    T1=np.array([fitT1(Lobs[q].mean(0)) for q in range(Q)])
    base=np.exp(-T[None,:]/T1[:,None])
    ind,ok=below(Lobs); cnt=ok.sum(1); ob=(ind&ok).sum(1)
    rec=dict(chip=run['chip'],run=run['run'],Q=Q,obs=ob.sum()/cnt.sum(),idle0=idle0,
             T1med=float(np.median(T1)))
    # published procedure on real data (reproduce)
    rec['pub_null']=published_null(Lobs,pr)
    rec['var']={}
    for name,v in V.items():
        b=idle0 if v['floor'] else 0.
        pop=base+b*(1-base)
        n=calib_n(Lobs,v['sep'],v['fd'],v['epp'],v['e0'],pop)
        Ls=chain(pop,n,v['sep'],v['fd'],v['epp'],v['e0'],R,v['shared'])
        bi,bo=below(Ls)
        p0=(bi&bo).sum((0,2))/np.maximum(bo.sum((0,2)),1)
        exp_ob=p0*cnt
        null=exp_ob.sum()/cnt.sum()
        idx=rng.integers(0,Q,(B,Q))
        o_b=ob[idx].sum(1)/cnt[idx].sum(1); n_b=exp_ob[idx].sum(1)/cnt[idx].sum(1)
        d=o_b-n_b
        lo,hi=np.percentile(d,[2.5,97.5]); olo,ohi=np.percentile(o_b,[2.5,97.5]); nlo,nhi=np.percentile(n_b,[2.5,97.5])
        rec['var'][name]=dict(n=n,null=null,d=rec['obs']-null,dlo=lo,dhi=hi,olo=olo,ohi=ohi,nlo=nlo,nhi=nhi,
                              verdict='ABOVE' if lo>0 else 'BELOW' if hi<0 else 'CONSISTENT')
        if name=='main':
            # calibration of published procedure on simulated pure-exponential data through chain
            cal=[published_null(Ls[i],pr) for i in range(6)]
            bi0,bo0=below(Ls[:6]); rec['sim_pub_proc']=float(np.mean(cal)); rec['sim_true_share']=float((bi0&bo0).sum()/bo0.sum())
    rec['final']=(lambda vs: vs[0] if len(set(vs))==1 else 'CANNOT_TELL')([x['verdict'] for x in rec['var'].values()])
    out.append(rec); print(rec['chip'],rec['run'],'obs %.3f pubnull %.3f main %.3f d[%.3f,%.3f] n=%d simpubproc %.3f simtrue %.3f %s'%(rec['obs'],rec['pub_null'],rec['var']['main']['null'],rec['var']['main']['dlo'],rec['var']['main']['dhi'],rec['var']['main']['n'],rec['sim_pub_proc'],rec['sim_true_share'],rec['final']),flush=True)
json.dump(out,open('results.json','w'),default=float)
