"""Outside referee check of the idle-in-0 offset (claim: 300 us level above the zero-wait 0 level; 600 us; per-job widths) and of the
decay-time repeats (shot-noise simulation of the median fading time). Written 22 Oct 2026 by an outside referee after the data. Run from this folder."""
import json,math,statistics as st,numpy as np
B='/Users/heuremen/bones-forge/'
d=json.load(open(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)),'data','time_runs_ruler_copies.json')))['runs']
print(len(d))
pos=0;n600up=0
for r in d:
    R=np.array(r['rows']);i3=R[:,13];i6=R[:,14]
    m3=100*i3.mean();m6=100*i6.mean()
    se3=100*i3.std(ddof=1)/math.sqrt(len(i3));
    # between-qubit scatter includes true per-qubit differences + shot noise
    print(r['chip'][4:9],r['run'],len(R),'idle300 %.2f se %.2f | idle600 %.2f diff %.2f | n>0 %.0f%%'%(m3,se3,m6,m6-m3,100*(i3>0).mean()))
    pos+=m3>0;n600up+=m6>m3
print('pos',pos,'600>300',n600up)
# tau noise
def tau(l): return -300/math.log(min(max(l,.02),.98))
import itertools
for chip in ['ibm_marrakesh','ibm_kingston','ibm_pittsburgh']:
  for r in [x for x in d if x['chip']==chip]:
    R=np.array(r['rows']); L=R[:,[3,7,11]]  # 300us levels a,b,c
    T=np.vectorize(tau)(L)
    med=np.median(T,axis=0)
    # per-qubit copy-to-copy scatter of level (within-qubit, shot noise + drift)
    sd=L.std(axis=1,ddof=1); 
    # level sd per copy estimated from within-qubit scatter
    s=math.sqrt((sd**2).mean())
    # bootstrap-free: simulate median change from independent noise s on levels
    rng=np.random.default_rng(0); mean=L.mean(1)
    rng_meds=[]
    for _ in range(400):
        Ln=mean+rng.normal(0,s,len(mean)); rng_meds.append(np.median(np.vectorize(tau)(Ln)))
    # copy-range: simulate three independent copies
    rngs=[]
    for _ in range(400):
        m=[np.median(np.vectorize(tau)(mean+rng.normal(0,s,len(mean)))) for c in range(3)]
        rngs.append(max(m)-min(m))
    print(chip[4:9],r['run'],len(R),'meds',np.round(med,1),'range %.1f'%(med.max()-med.min()),'| within-qubit level sd %.4f; sim median-sd %.2f, sim range of 3 mean %.1f (p95 %.1f)'%(s,np.std(rng_meds),np.mean(rngs),np.percentile(rngs,95)),'mean L300 %.3f'%L.mean())
z=[];dd=[];ch={}
for r in d:
    R=np.array(r['rows']);a=R[:,13];b=R[:,14]
    z.append(a.mean()/(a.std(ddof=1)/math.sqrt(len(a))));dd.append(100*(b-a).mean());ch.setdefault(r['chip'],[]).append(100*(b-a).mean())
print('z',np.round(z,1),'n z<2',sum(x<2 for x in z),'<3',sum(x<3 for x in z));print('mean 600-300 %.2f'%np.mean(dd));print({k:round(float(np.mean(v)),2) for k,v in ch.items()})
