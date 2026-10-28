from sim import *
for run in [r for r in runs if (r['chip'],r['run']) in [('ibm_marrakesh',301),('ibm_pittsburgh',361),('ibm_aachen',871),('ibm_fez',401)]]:
    A=np.array(run['rows']);Q=len(A)
    Lobs=np.stack([A[:,1:5],A[:,5:9],A[:,9:13]],1)
    T1=np.array([fitT1(Lobs[q].mean(0)) for q in range(Q)]); pop=np.exp(-T[None,:]/T1[:,None])
    row=[]
    for sep in (4,6,9):
        Ls=chain(pop,300,sep,.03,.01,.005,200,mode='bits'); bi,bo=below(Ls); row.append('sep%d raw-bit share %.1f%%'%(sep,100*(bi&bo).sum()/bo.sum()))
        Ls=chain(pop,300,sep,.03,.01,.005,200,mode='iq'); bi,bo=below(Ls); row.append('(IQ-level %.1f%%)'%(100*(bi&bo).sum()/bo.sum()))
    print(run['chip'],run['run'],' '.join(row))
