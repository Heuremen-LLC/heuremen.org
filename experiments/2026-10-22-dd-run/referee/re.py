import numpy as np, json, sys
chip=sys.argv[1]; D='/Users/heuremen/bones-forge/arxiv_one_readout/dd_run/'
z=np.load(D+'raw_dd_%s.npz'%chip); j=json.load(open(D+'job_dd_%s.json'%chip)); meta=j['meta']; S=np.array(j['info']['S']); n=j['n_qubits']
def ex(dd,mode,N,st,copies=(0,1)):
    acc=[]
    for i,m in enumerate(meta):
        if m['dd']==dd and m['mode']==mode and m['N']==N and m['setting']==st and m['copy'] in copies:
            b=np.unpackbits(z['fin_%d'%i],axis=1,bitorder='big')[:,::-1][:,:n]; acc.append(1-2*b.astype(int))
    return np.concatenate(acc).mean(0)[S]
ph=lambda d,mo,N: np.arctan2(ex(d,mo,N,'YY'),ex(d,mo,N,'XX'))
co=lambda d,mo,N: np.hypot(ex(d,mo,N,'YY'),ex(d,mo,N,'XX'))
wrap=lambda x:(x+np.pi)%(2*np.pi)-np.pi
Ns=[2,4,8,16]
def gap(d):
    g={N:wrap(ph(d,'M',N)-ph(d,'D',N)) for N in Ns}; t={2:g[2]}
    for a,b in zip(Ns[:-1],Ns[1:]): t[b]=t[a]+wrap(g[b]-g[a])
    return np.degrees(np.stack([t[N] for N in Ns],1))
X=np.stack([np.array(Ns,float),np.ones(4)],1)
sl=lambda d: np.linalg.lstsq(X,gap(d).T,rcond=None)[0][0]
fl=5/np.sqrt(8000)
alive=lambda d: np.all([co(d,mo,N)>=fl for N in Ns for mo in 'DM'],0)
lo,on=alive('off'),alive('on'); both=lo&on
print(chip,'R1',np.median(sl('off')[lo]),'R2',np.median(sl('on')[on]),'R3',np.median((sl('on')-sl('off'))[both]))
for N in Ns:
    r4=co('off','D',N)-co('off','M',N); r5=co('on','D',N)-co('on','M',N)
    print(' N',N,'R4 %+.4f R5 %+.4f R6 %+.4f R7 %+.4f'%(np.median(r4),np.median(r5),np.median(r5-r4),np.median(co('on','D',N)-co('off','D',N))))
# mean-based check, and raw phase of D arms (is the idle-arm phase itself changed by echo?)
print(' median phase D off/on N=8 (deg):',np.degrees(np.median(ph('off','D',8))),np.degrees(np.median(ph('on','D',8))))
print(' mean coh D off N=0, N=8; on N=0,N=8: %.3f %.3f %.3f %.3f'%(co('off','D',0).mean(),co('off','D',8).mean(),co('on','D',8).mean(),co('on','D',16).mean()))
print(' mean coh M off N=8 %.3f on %.3f ; N=16 D off %.3f M off %.3f D on %.3f M on %.3f'%(co('off','M',8).mean(),co('on','M',8).mean(),co('off','D',16).mean(),co('off','M',16).mean(),co('on','D',16).mean(),co('on','M',16).mean()))
# mean-of-qubits R6 (not median) with qubit bootstrap
rng=np.random.default_rng(1)
for N in (2,4,8):
    d=(co('on','D',N)-co('on','M',N))-(co('off','D',N)-co('off','M',N)); bm=[rng.choice(d,len(d)).mean() for _ in range(2000)]
    r4=co('off','D',N)-co('off','M',N)
    print(' N',N,'R6 mean %+.4f CI [%+.4f,%+.4f]; R4 median %+.4f -> halfwidth/R4 = %.2f'%(d.mean(),*np.percentile(bm,[2.5,97.5]),np.median(r4),(np.percentile(bm,97.5)-np.percentile(bm,2.5))/2/np.median(r4)))
