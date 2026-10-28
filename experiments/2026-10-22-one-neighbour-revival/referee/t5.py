import sys;sys.path.insert(0,'/private/tmp/claude-501/-Users-heuremen/80087010-0cb3-4ae5-9e13-362e6d5086ad/scratchpad/ref_zz')
from lib import *
G=np.arange(0.10,2.5001,0.0025)
for chip in('ibm_phoenix','ibm_fez'):
    Jb=J(chip);taus=np.array(Jb.taus);ts=taus*1e-6
    C={a:np.array([coh(Jb,a,t,(0,1)) for t in Jb.taus]) for a in('Z','P1','P2')}
    rows=[]
    for j,(b1,a,b2) in enumerate(Jb.tr):
        for arm,q in(('P1',b1),('P2',b2)):
            f=Jb.zz[tuple(sorted((a,q)))]
            if f<=0:continue
            mod=C['Z'][:,j][None]*np.abs(np.cos(np.pi*G[:,None]*f*ts[None]));k=np.argmin(((mod-C[arm][:,j][None])**2).sum(1))
            jq=next(g['value'] for g in Jb.p['general'] if g['name']=='jq_%d%d'%(min(a,q),max(a,q)) or g['name']=='jq_%d%d'%(max(a,q),min(a,q)))
            rows.append((f/1e3,G[k],G[k]*f/1e3,jq*1e6,1e6*0 ))
    r=np.array(rows);print(chip,'pairs',len(r))
    print(' corr(zz, fitted beat s*zz)',np.corrcoef(r[:,0],r[:,2])[0,1],' spearman-ish',np.corrcoef(np.argsort(np.argsort(r[:,0])),np.argsort(np.argsort(r[:,2])))[0,1])
    for lo,hi in((0,5),(5,10),(10,20),(20,30),(30,60)):
        m=(r[:,0]>=lo)&(r[:,0]<hi)
        if m.sum():print('  zz %d-%d kHz n=%d median s=%.2f  [per-pair s: %s]'%(lo,hi,m.sum(),np.median(r[m,1]),np.round(np.sort(r[m,1]),2)[:12]))
    print(' jq (kHz) median',np.median(r[:,3]),'max',r[:,3].max())
    # zz/jq^2 ratio
    print(' zz vs jq^2 corr',np.corrcoef(r[:,0],r[:,3]**2)[0,1])
