import sys;sys.path.insert(0,'/Users/heuremen/bones-forge/arxiv_one_readout/zz_revival')
import analyze_zz as A, numpy as np
def run(folder):
    A.HERE=folder
    Jb=A.Job('ibm_fez');T=np.array(Jb.taus);m=Jb.m
    C={a:np.array([A.coh(Jb,a,t,(0,1))[0] for t in Jb.taus]) for a in('Z','P1')}
    R,info=A.readings(Jb,(0,1))
    h=R['R2 return height at revival / Z (as measured)'];d=R['R2 return height at revival / Z, minus ZZ+T1 prediction (1+exp(-t*/T1_B1))/2']
    s=R['R1 fitted scale s (P1, B1 alone in +)'];f1=info['f1'];ts=1/(s*f1)*1e6
    zc=np.full(m,np.nan)
    for j in range(m):
        if np.isfinite(h[j]):k=np.argsort(abs(T-ts[j]))[:3];zc[j]=C['Z'][k,j].mean()
    ok=np.isfinite(d)
    r=lambda mask:np.nanmedian(np.where(mask&ok,d,np.nan))
    return ok.sum(),(zc[ok]<=0.15).sum(),r(np.ones(m,bool)),r(zc>0.15),r(zc<=0.15),np.nanmedian(np.where((zc<=0.15)&ok,h,np.nan)),np.nanmedian(zc[ok])
print('folder n  n(Z<=.15)  R2all  R2(Z>.15)  R2(Z<=.15)  h(Z<=.15)  median Z(t*)')
print('REAL',*np.round(run('/Users/heuremen/bones-forge/arxiv_one_readout/zz_revival'),3))
for d in['sg_11','sg_12','sg_13','sg_14','sg_15','sg_16','sg_18','sg_19','sg_20','sf_1.0_1','sf_1.0_2']:
    print(d,*np.round(run('/private/tmp/claude-501/-Users-heuremen/80087010-0cb3-4ae5-9e13-362e6d5086ad/scratchpad/ref_zz/'+d),3))
