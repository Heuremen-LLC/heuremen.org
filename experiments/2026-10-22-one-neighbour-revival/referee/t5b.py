import sys;sys.path.insert(0,'/private/tmp/claude-501/-Users-heuremen/80087010-0cb3-4ae5-9e13-362e6d5086ad/scratchpad/ref_zz')
from lib import *
Jb=J('ibm_phoenix');taus=Jb.taus
C={a:np.array([coh(Jb,a,t,(0,1)) for t in taus]) for a in('Z','P1','P2')}
np.set_printoptions(linewidth=250)
print('taus',taus[:41])
for j in (17,3,2,5):
    b1,a,b2=Jb.tr[j];print('triple',Jb.tr[j],'zz1 %.1f zz2 %.1f kHz'%(Jb.zz[tuple(sorted((a,b1)))]/1e3,Jb.zz[tuple(sorted((a,b2)))]/1e3))
    print(' Z ',np.round(C['Z'][:,j],2));print(' P1/Z',np.round(C['P1'][:,j]/C['Z'][:,j],2));print(' P2/Z',np.round(C['P2'][:,j]/C['Z'][:,j],2))
