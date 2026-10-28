import sys;sys.path.insert(0,'/private/tmp/claude-501/-Users-heuremen/80087010-0cb3-4ae5-9e13-362e6d5086ad/scratchpad/ref_zz')
from lib import *
for chip in('ibm_phoenix','ibm_fez'):
    Jb=J(chip);ts=np.array(Jb.taus)*1e-6
    C={a:np.array([coh(Jb,a,t,(0,1)) for t in Jb.taus]) for a in('Z','P1','P2')}
    F=np.arange(1e3,100e3,250.);fb=[]
    for j in range(Jb.m):
        r=[]
        for arm in('P1','P2'):
            mod=C['Z'][:,j][None]*np.abs(np.cos(np.pi*F[:,None]*ts[None]));r.append(F[np.argmin(((mod-C[arm][:,j][None])**2).sum(1))]/1e3)
        fb.append(r)
    fb=np.array(fb);z=np.array([[Jb.zz[tuple(sorted((a,b1)))]/1e3,Jb.zz[tuple(sorted((a,b2)))]/1e3] for b1,a,b2 in Jb.tr])
    ok=(fb<99).all(1)&(z>0).all(1)
    print(chip,'P1 vs P2 free-fit beat (kHz) corr %.2f; props zz1 vs zz2 corr %.2f'%(np.corrcoef(fb[ok,0],fb[ok,1])[0,1],np.corrcoef(z[ok,0],z[ok,1])[0,1]))
    print(' within-triple |log(f1/f2)| median fitted %.2f vs props %.2f'%(np.median(abs(np.log(fb[ok,0]/fb[ok,1]))),np.median(abs(np.log(z[ok,0]/z[ok,1])))))
