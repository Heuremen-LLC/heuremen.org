import json,sys
sys.path.insert(0,'.')
for chip in ('ibm_phoenix','ibm_fez'):
    rec=json.load(open(f'job_zz_{chip}.json')); props=json.load(open(f'props_zz_{chip}.json'))
    edges={tuple(sorted(e)) for e in rec['edges']}
    names=[g for g in props['general'] if g['name'].startswith('zz_')]
    print(chip,'edges',len(edges),'zz names',len(names), 'sample', names[:2])
    cnt={}
    diff=[]
    for g in names:
        d=g['name'][3:]
        naive=[(int(d[:k]),int(d[k:])) for k in range(1,len(d)) if (min(int(d[:k]),int(d[k:])),max(int(d[:k]),int(d[k:]))) in edges]
        fixed=sorted({(min(int(d[:k]),int(d[k:])),max(int(d[:k]),int(d[k:]))) for k in range(1,len(d)) if all(x=='0' or x[0]!='0' for x in (d[:k],d[k:])) and (min(int(d[:k]),int(d[k:])),max(int(d[:k]),int(d[k:]))) in edges})
        nv={(min(a,b),max(a,b)) for a,b in naive}
        if len(naive)!=1 or len(fixed)!=1 or nv!=set(fixed): diff.append((g['name'],naive,fixed,g['value']))
        for e in fixed: cnt.setdefault(e,[]).append((g['name'],g['value']))
    print(' differing names:'); [print('  ',x) for x in diff]
    print(' edges with no zz:',[e for e in edges if e not in cnt]); print(' edges with >1 zz entries:',{e:v for e,v in cnt.items() if len(v)>1})
    print(' zz entries mapping to non-edge/none:',[g['name'] for g in names if not any(g['name']==n for e in cnt for n,_ in cnt[e])])
    for t in rec['info']['triples']:
        b1,a,b2=t
        for q in (b1,b2):
            e=tuple(sorted((a,q)))
            if e not in edges or e not in cnt: print(' TRIPLE PROBLEM',t,e)
    # names for all triple edges
    print(' triple values kHz:',[(t, [round(abs(cnt[tuple(sorted((t[1],q)))][0][1])*1e6,2) for q in (t[0],t[2])]) for t in rec['info']['triples']][:40])
