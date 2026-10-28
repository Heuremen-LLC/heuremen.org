"""Fez-like synthetic: real fez triples, real props (zz, T1, T2, readout), fez tau grid; same physics as synth_zz.py. python synth_fez.py OUTDIR SCALE SEED"""
import json,os,sys,shutil
import numpy as np
D='/Users/heuremen/bones-forge/arxiv_one_readout/zz_revival'; sys.path.insert(0,'/private/tmp/claude-501/-Users-heuremen/80087010-0cb3-4ae5-9e13-362e6d5086ad/scratchpad/ref_zz')
from lib import zzmap,qp
OUT,S,SEED=sys.argv[1],float(sys.argv[2]),int(sys.argv[3]);os.makedirs(OUT,exist_ok=True);rng=np.random.default_rng(SEED)
rec=json.load(open(D+'/job_zz_ibm_fez.json'));props=json.load(open(D+'/props_zz_ibm_fez.json'))
edges={tuple(sorted(e)) for e in rec['edges']};zzr=zzmap(props,edges)
trip=[tuple(t) for t in rec['info']['triples']];m=len(trip);n=rec['n_qubits'];SHOTS=1000
TAUS=[3+7.5*k for k in range(41)];ARMS=('Z','P1','P2','P12','O1','EZ','EP1')
T1=np.array([qp(props,q,'T1') for q in range(n)]);T2=np.array([qp(props,q,'T2') for q in range(n)]);ro=np.array([qp(props,q,'readout_error') for q in range(n)])
det=rng.uniform(-20e3,20e3,n);zz=zzr
shutil.copy(D+'/props_zz_ibm_fez.json',OUT+'/props_zz_ibm_fez.json')
def pack(b):
    nb=(b.shape[1]+7)//8;pad=np.zeros((b.shape[0],nb*8),np.uint8);pad[:,:b.shape[1]]=b;return np.packbits(pad[:,::-1],axis=1,bitorder='big')
meta=[];out={};i=0
for c in(0,1):
  for arm in ARMS:
    for tau in TAUS:
      for st in 'XY':
        t=tau*1e-6;A=np.zeros((SHOTS,m),np.uint8);b1r=np.zeros_like(A);b2r=np.zeros_like(A)
        for j,(b1,a,b2) in enumerate(trip):
            echo=arm in('EZ','EP1')
            def branch(on,Tq):
                st1=(rng.random(SHOTS)<0.5) if on=='+' else (np.ones(SHOTS,bool) if on=='1' else np.zeros(SHOTS,bool))
                td=rng.exponential(Tq*1e-6,SHOTS);t_in1=np.where(st1,np.minimum(td,t),0.0)
                if echo:
                    h=t/2;first=np.where(st1,np.minimum(td,h),0.0);t_in1=first-(t_in1-first)
                return t_in1,st1&(td>t)
            s1='+' if arm in('P1','P12','EP1') else '1' if arm=='O1' else '0';s2='+' if arm in('P2','P12') else '0'
            ti1,e1=branch(s1,T1[b1]);ti2,e2=branch(s2,T1[b2])
            k1=zz[tuple(sorted((a,b1)))];k2=zz[tuple(sorted((a,b2)))]
            phi=2*np.pi*(S*k1*ti1+S*k2*ti2)+(0 if echo else 2*np.pi*det[a]*t)+rng.normal(0,1,SHOTS)*np.sqrt(2*t/(T2[a]*1e-6))
            phi=phi+(0 if st=='X' else -np.pi/2)
            A[:,j]=(rng.random(SHOTS)<(1-np.cos(phi))/2)
            b1r[:,j]=e1^(rng.random(SHOTS)<ro[b1]);b2r[:,j]=e2^(rng.random(SHOTS)<ro[b2])
        out['a_%d'%i]=pack(A);out['b1_%d'%i]=pack(b1r);out['b2_%d'%i]=pack(b2r)
        meta.append(dict(arm=arm,tau_us=tau,tau_dt=0,setting=st,copy=c));i+=1
np.savez_compressed(OUT+'/raw_zz_ibm_fez.npz',**out)
json.dump(dict(run='zz',chip='ibm_fez',n_qubits=n,job_id='SYNTHFEZ-%s-%d'%(S,SEED),meta=meta,shots=SHOTS,info=dict(triples=trip),edges=sorted(edges)),open(OUT+'/job_zz_ibm_fez.json','w'))
