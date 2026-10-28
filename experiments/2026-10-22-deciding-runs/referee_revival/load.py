import json,numpy as np
D='/Users/heuremen/bones-forge/arxiv_one_readout/deciding_runs/'
def bits(arr,n): return np.unpackbits(np.asarray(arr,np.uint8),axis=1,bitorder='big')[:,::-1][:,:n].astype(np.int8)
def load(chip):
    rec=json.load(open(D+f'job_revival_ibm_{chip}.json')); z=np.load(D+f'raw_revival_ibm_{chip}.npz'); props=json.load(open(D+f'props_revival_ibm_{chip}.json'))
    A=rec['info']['A']; B=rec['info']['B']; meta=rec['meta']; taus=np.array(sorted({m['tau'] for m in meta}),float)
    # per (copy,arm,tau,setting) -> mean_pm array over A
    M={}
    for i,m in enumerate(meta):
        M[(m['copy'],m['arm'],m['tau'],m['setting'])]=np.mean(1-2*bits(z['a_%d'%i],len(A)),0)
    return rec,props,A,B,taus,M
def curves(M,taus,arm,copies=(0,1)):
    # returns X,Y,C arrays [nA, ntau] pooled over copies (equal shots)
    X=np.stack([np.mean([M[(c,arm,int(t),'X')] for c in copies],0) for t in taus],1)
    Y=np.stack([np.mean([M[(c,arm,int(t),'Y')] for c in copies],0) for t in taus],1)
    return X,Y,np.hypot(X,Y)
def edges(props):
    s=set()
    for g in props['gates']:
        if g['gate']=='cz': s.add(tuple(sorted(g['qubits'])))
    return s
