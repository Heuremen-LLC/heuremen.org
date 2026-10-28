import json,numpy as np,os
D='/Users/heuremen/bones-forge/arxiv_one_readout/zz_revival'
class J:
    def __init__(s,chip):
        s.r=json.load(open(f'{D}/job_zz_{chip}.json')); s.z=np.load(f'{D}/raw_zz_{chip}.npz'); s.p=json.load(open(f'{D}/props_zz_{chip}.json'))
        s.tr=[tuple(t) for t in s.r['info']['triples']]; s.m=len(s.tr); s.meta=s.r['meta']
        s.taus=sorted({m['tau_us'] for m in s.meta}); s.edges={tuple(sorted(e)) for e in s.r['edges']}
        s.ix={}
        for i,m in enumerate(s.meta): s.ix.setdefault((m['arm'],m['tau_us'],m['setting'],m['copy']),i)
        s.zz=zzmap(s.p,s.edges)
    def bit(s,reg,i):
        a=s.z[f'{reg}_{i}']; nb=a.shape[1]
        return np.stack([(a[:,nb-1-j//8]>>(j%8))&1 for j in range(s.m)],1)   # independent unpack: qubit j = LSB of last byte
def zzmap(p,edges):
    out={}
    for g in p['general']:
        if not g['name'].startswith('zz_'): continue
        d=g['name'][3:]; c=[]
        for k in range(1,len(d)):
            a,b=d[:k],d[k:]
            if (len(a)>1 and a[0]=='0') or (len(b)>1 and b[0]=='0'): continue
            e=tuple(sorted((int(a),int(b))))
            if e in edges and e not in c: c.append(e)
        assert len(c)==1,(g['name'],c); out[c[0]]=abs(g['value'])*1e9
    return out
def qp(p,q,n):
    return next((x['value'] for x in p['qubits'][q] if x['name']==n),(0.02 if n=='readout_error' else 100.0))
def expv(J_,arm,tau,st,copies,cond=None):
    num=np.zeros(J_.m);den=np.zeros(J_.m)
    for c in copies:
        i=J_.ix[(arm,tau,st,c)]; v=1-2*J_.bit('a',i).astype(float)
        w=np.ones_like(v) if cond is None else (J_.bit(cond[0],i)==cond[1]).astype(float)
        num+=(v*w).sum(0);den+=w.sum(0)
    return num/np.maximum(den,1),den
def coh(J_,arm,tau,copies):
    x,_=expv(J_,arm,tau,'X',copies);y,_=expv(J_,arm,tau,'Y',copies);return np.hypot(x,y)
