"""Hand-added counts in the paper/site, recomputed from published files. stdlib + numpy only. Read-only."""
import json, re, glob, os, collections, datetime as dt
import numpy as np
B = '/Users/heuremen/bones-forge'
X = B + '/heuremen.org/experiments'
A = B + '/arxiv_one_readout'
def verdict(ok, mine, paper):
    print('   -> %s   (computed: %s | paper: %s)' % ('REPRODUCES' if ok else 'DIFFERS', mine, paper))
def hdr(s): print('\n' + '=' * 100 + '\n' + s)

ibm = json.load(open(X + '/job-index/ibm_jobs_2026-09-24_to_10-22.json'))
rig = json.load(open(X + '/job-index/rigetti_tasks_2026-09-25_to_10-05.json'))
ix = {j['job_id']: j for j in ibm}
day = lambda j: j['created_et'][:10]
byday = lambda d: sorted([j for j in ibm if day(j) == d], key=lambda j: j['created_et'])
bk = lambda js: dict(sorted(collections.Counter(j['backend'].replace('ibm_', '') for j in js).items()))

hdr('1a. "25--28 Sep: Rigetti, 563 tasks"  (US Eastern day = created UTC - 4 h, as job-index/count_jobs.py)')
et = lambda s: (dt.datetime.fromisoformat(s) - dt.timedelta(hours=4)).date().isoformat()
rd = collections.Counter(et(t['created']) for t in rig)
sel = {d: rd[d] for d in ('2026-09-25', '2026-09-26', '2026-09-28')}
print('   Rigetti tasks per Eastern day (all days in file):', dict(sorted(rd.items())))
print('   25, 26, 28 Sep:', sel)
verdict(sum(sel.values()) == 563, sum(sel.values()), 563)

hdr('1b. "1--2 Oct: kingston, phoenix, pittsburgh, boston (12)" = all 1 Oct jobs + first 9 of 2 Oct by created_et')
o1, o2 = byday('2026-10-01'), byday('2026-10-02')
first12 = o1 + o2[:9]
print('   1 Oct jobs: %d %s' % (len(o1), bk(o1)))
print('   first 9 of 2 Oct: %s' % bk(o2[:9]))
print('   10th 2 Oct job backend:', o2[9]['backend'])
print('   the 12 together: %d %s' % (len(first12), bk(first12)))
verdict(len(first12) == 12 and set(bk(first12)) == {'kingston', 'phoenix', 'pittsburgh', 'boston'},
        '%d %s' % (len(first12), list(bk(first12))), '12; kingston, phoenix, pittsburgh, boston')

hdr('1c. "2--7 Oct: 7 processors (305): pittsburgh 86, kingston 75, marrakesh 55, aachen 25, boston 24, fez 20, phoenix 20; 63 on 2 Oct, 43 on 4, 26 on 5, 48 on 6, 125 on 7; three errored"')
firstids = {j['job_id'] for j in first12}
r = [j for j in ibm if '2026-10-02' <= day(j) <= '2026-10-07' and j['job_id'] not in firstids]
pc = collections.Counter(j['backend'].replace('ibm_', '') for j in r)
dc = collections.Counter(day(j) for j in r)
paper_p = dict(pittsburgh=86, kingston=75, marrakesh=55, aachen=25, boston=24, fez=20, phoenix=20)
paper_d = {'2026-10-02': 63, '2026-10-04': 43, '2026-10-05': 26, '2026-10-06': 48, '2026-10-07': 125}
err = [j for j in r if j['status'] != 'DONE']
print('   total:', len(r), '| processors:', dict(pc.most_common()))
print('   per Eastern day:', dict(sorted(dc.items())))
print('   status != DONE:', [(day(j), j['backend'], j['status']) for j in err])
verdict(len(r) == 305, len(r), 305)
verdict(dict(pc) == paper_p, dict(pc.most_common()), paper_p)
verdict(dict(dc) == paper_d, dict(sorted(dc.items())), paper_d)
verdict(len(err) == 3 and all(j['status'] == 'ERROR' for j in err), '%d non-DONE: %s' % (len(err), [(day(j)[5:], j['backend'].replace('ibm_', ''), j['status']) for j in err]), 'three errored, one each on 2, 4, 6 Oct')
allerr = [(day(j), j['backend'], j['status']) for j in ibm if j['status'] != 'DONE']
print('   (all non-DONE jobs in the whole index, for context: %s)' % allerr)

hdr('1d. "7 Oct (within the 305): pittsburgh 9, fez 3, aachen 2, kingston 18 (32)"  jobs named in 2026-10-07-equal-idle/results_oct07.text')
txt = open(X + '/2026-10-07-equal-idle/results_oct07.text').read()
ids = []
for line in txt.splitlines():
    m = re.search(r'\| jobs (.*)$', line)
    if m:
        for i in m.group(1).split(','): ids.append(i.strip())
u = list(dict.fromkeys(ids))
print('   job ids listed (with repeats across families that share jobs): %d; distinct: %d' % (len(ids), len(u)))
miss = [i for i in u if i not in ix]
print('   not in job index:', miss)
sub = [ix[i] for i in u if i in ix]
print('   days of these jobs:', dict(collections.Counter(day(j) for j in sub)))
print('   by backend:', bk(sub))
verdict(bk(sub) == dict(aachen=2, fez=3, kingston=18, pittsburgh=9) and len(sub) == 32, '%d %s' % (len(sub), bk(sub)), '32: pittsburgh 9, fez 3, aachen 2, kingston 18')
print('   (README.md.text of that folder also states the same 32 / 18 / 9 / 3 / 2 in words; the count above is from the job ids.)')
inrest = [j for j in sub if j['job_id'] in {x['job_id'] for x in r}]
print('   all %d are among the 305 of 1c: %s' % (len(sub), len(inrest) == len(sub)))

hdr('1e. "24 Sep: fez 3, kingston 4, marrakesh 41 (48; 17 in the run log)"')
s24 = byday('2026-09-24')
print('   job index 24 Sep:', len(s24), bk(s24))
verdict(len(s24) == 48 and bk(s24) == dict(fez=3, kingston=4, marrakesh=41), '%d %s' % (len(s24), bk(s24)), '48: fez 3, kingston 4, marrakesh 41')
found = {}
for f in glob.glob(X + '/2026-09-24-temporal-spatial-correlations/**/*', recursive=True):
    if os.path.isfile(f):
        try: t = open(f).read()
        except Exception: continue
        for i in re.findall(r'\b[a-z0-9]{20}\b', t):
            if i in ix: found.setdefault(i, set()).add(os.path.basename(f))
print('   NO FILE IS A "RUN LOG" of 24 Sep: grep of experiments/ for run log / run_log finds none for that day.')
print('   Plain count of distinct 24 Sep job ids that the published 24 Sep folder names: %d %s' % (len(found), bk([ix[i] for i in found])))
print('   -> NO DATA for "17 in the run log" (the log is not published). The %d above is a different quantity (ids named in saved result files), not a substitute; it does not equal 17.' % len(found))

hdr('1f. "21 Oct: fez, kingston, marrakesh, pittsburgh (4; a boston job errored)"  jobs in 2026-10-21-readout-vs-idle/job_*.json')
js = []
for f in sorted(glob.glob(X + '/2026-10-21-readout-vs-idle/job_*.json')):
    t = open(f).read()
    for i in set(re.findall(r'\b[a-z0-9]{20}\b', t)):
        if i in ix: js.append(ix[i])
print('   jobs:', [(j['backend'].replace('ibm_', ''), j['status'], j['created_et'][11:19]) for j in js])
ok = sorted(j['backend'] for j in js if j['status'] == 'DONE'); bad = [j['backend'] for j in js if j['status'] != 'DONE']
verdict(ok == ['ibm_fez', 'ibm_kingston', 'ibm_marrakesh', 'ibm_pittsburgh'] and bad == ['ibm_boston'], 'DONE %s; not done %s' % (ok, bad), '4 DONE (fez, kingston, marrakesh, pittsburgh); boston errored')

# ---------------------------------------------------------------- item 2
hdr('2. 7 Oct equal-idle arms (results_oct07.json). Arm = (processor, family, pair) with a "kick all" reading, excluding POOL lines and arms with n=1 (CANNOT TELL, no median over qubits).')
R = json.load(open(X + '/2026-10-07-equal-idle/results_oct07.json'))
arms = collections.OrderedDict()
for c, f, l, k, p in R: arms.setdefault((c, f, l), {})[k] = p
tab = {a: v for a, v in arms.items() if a[1] != 'POOL' and 'kick all' in v and v['kick all'] is not None and v['kick all'][3] > 1}
print('   arms in table:', len(tab), '(paper: "25 arms"; the table tab:oct07 has 25 rows)')
out = lambda v, k: (v[k][4] if v.get(k) is not None else 'COPIES SHIFTED')
k25 = {a: v for a, v in tab.items() if a[0] == 'ibm_kingston' and 'CLOCK' not in a[1]}
print('   kingston 25-measured arms (state-dark 3, spaced 3, gapped 3, wire paint 3, carriers 2):', len(k25))
b_all = sum(out(v, 'kick all') == 'BELOW' for v in k25.values())
b_d2 = sum(out(v, 'kick d<=2') == 'BELOW' for v in k25.values())
verdict(len(k25) == 14 and b_all == 13, '%d arms, kick-all BELOW in %d' % (len(k25), b_all), '14 arms, below zero in 13')
medall = [v['kick all'][0] for v in k25.values() if out(v, 'kick all') == 'BELOW']
print('   medians of those BELOW arms: %.2f to %.2f' % (min(medall), max(medall)))
verdict(b_d2 == 12, 'kick d<=2 BELOW in %d of %d' % (b_d2, len(k25)), '12 of 14')
d5 = collections.Counter(out(v, 'kick d>=5') for v in tab.values())
print('   kick d>=5 outcomes over the %d table arms: %s (paper: none BELOW; reversed-order arm the only ABOVE)' % (len(tab), dict(d5)))
verdict(d5.get('BELOW', 0) == 0 and d5.get('ABOVE', 0) == 1, 'BELOW %d, ABOVE %d' % (d5.get('BELOW', 0), d5.get('ABOVE', 0)), 'no BELOW; one ABOVE')
co = collections.Counter(out(v, 'coh all') for v in tab.values())
print('   "coherence slope over all untouched qubits" read as the table column "coh all" outcome (BELOW / CANNOT TELL (incl. copies shifted) / CONSISTENT):', dict(co))
below = [a for a, v in tab.items() if out(v, 'coh all') == 'BELOW']
print('   BELOW arms:', [(a[0].replace('ibm_', ''), a[1], a[2]) for a in below])
ct = co.get('CANNOT TELL', 0) + co.get('COPIES SHIFTED', 0)
verdict((co.get('BELOW', 0), ct, co.get('CONSISTENT', 0)) == (3, 5, 17), 'below %d, cannot tell %d, consistent %d (above %d)' % (co.get('BELOW', 0), ct, co.get('CONSISTENT', 0), co.get('ABOVE', 0)), 'below 3, cannot tell 5, consistent 17 of 25')
n4 = {a: v['cohN4'] for a, v in tab.items() if v.get('cohN4') is not None and v['cohN4'][4] != 'CANNOT TELL' and v['cohN4'][0] == v['cohN4'][0]}
neg = sum(p[0] < 0 for p in n4.values())
print('   cohN4 (coherence difference at N=4): arms with a reading %d (the rest: copies shifted / cannot tell: %d); median < 0 in %d' % (len(n4), len(tab) - len(n4), neg))
print('   cohN4 BELOW arms:', [(a[0].replace('ibm_', ''), a[2]) for a, p in n4.items() if p[4] == 'BELOW'])
verdict((len(n4), neg) == (22, 20), '%d of %d negative' % (neg, len(n4)), '20 of the 22 arms with a reading')
print('   NOTE: "coherence slope ... below zero in 3 of 25" and "at N=4" are read as the table\'s "coh all" and the json\'s "cohN4" keys; the plain reading, not tuned.')

# ---------------------------------------------------------------- item 3
hdr('3. "echo and coherence-echo analyses read 20 of 396 (5.1%) and 32 of 528 (6.1%)"')
def tallies(path):
    tot = np.zeros(4, int)
    for l in open(path):
        m = re.match(r'seed \d+: TALLY readings (\d+): ABOVE (\d+), BELOW (\d+), CONSISTENT (\d+), CANNOT TELL (\d+)', l)
        if m: tot += np.array([int(x) for x in m.groups()[1:]])
    return tot
t1 = tallies(X + '/2026-10-22-dd-run/null_floor.text')
print('   echo run, null_floor.text (heuremen.org): readings %d, ABOVE %d, BELOW %d -> resolved %d (%.1f%%)' % (t1.sum(), t1[0], t1[1], t1[0] + t1[1], 100 * (t1[0] + t1[1]) / t1.sum()))
verdict((t1.sum(), t1[0] + t1[1]) == (396, 20), '%d of %d' % (t1[0] + t1[1], t1.sum()), '20 of 396')
p3 = A + '/dd_coherence/null_floor.txt'
t2 = tallies(p3)
t2b = np.zeros(4, int); n2b = 0
for f in sorted(glob.glob(A + '/dd_coherence/nulls/null_*.txt')):
    for l in open(f):
        m = re.match(r'TALLY readings (\d+): ABOVE (\d+), BELOW (\d+), CONSISTENT (\d+), CANNOT TELL (\d+)', l)
        if m: t2b += np.array([int(x) for x in m.groups()[1:]]); n2b += int(m.group(1))
print('   coherence-echo run (arxiv_one_readout/dd_coherence; NOT under heuremen.org/experiments):')
print('     null_floor.txt per-seed lines: readings %d, resolved %d (%.1f%%)' % (t2.sum(), t2[0] + t2[1], 100 * (t2[0] + t2[1]) / t2.sum()))
print('     recount from nulls/null_11..16.txt TALLY lines: readings %d, resolved %d (%.1f%%)' % (n2b, t2b[0] + t2b[1], 100 * (t2b[0] + t2b[1]) / n2b))
verdict((t2.sum(), t2[0] + t2[1], n2b, t2b[0] + t2b[1]) == (528, 32, 528, 32), '%d of %d (null_floor.txt), %d of %d (recount)' % (t2[0] + t2[1], t2.sum(), t2b[0] + t2b[1], n2b), '32 of 528')
print('   (counted from saved null outputs; the synthetic script was not re-run. The 528 file is in the arXiv folder, not in heuremen.org/experiments.)')

# ---------------------------------------------------------------- item 4
hdr('4. "the interval ends reach 0.5 to 1.3 times the cost itself" (echo run, coherence), results_dd.text [all]')
cur = None; Rv = {}
for l in open(X + '/2026-10-22-dd-run/results_dd.text'):
    m = re.match(r'dd (ibm_\w+) job', l)
    if m: cur = m.group(1); continue
    m = re.match(r'\s+(R4|R6) .*?N=\s*(\d+)(?: \(D - M\))? \[all\]\s+([+-][\d.e+-]+) \[([+-][\d.e+-]+), ([+-][\d.e+-]+)\]', l)
    if m and cur: Rv[(cur, m.group(1), int(m.group(2)))] = tuple(float(m.group(i)) for i in (3, 4, 5))
chips = ['ibm_fez', 'ibm_marrakesh', 'ibm_kingston']
paperN = {'ibm_fez': (2, 4, 8), 'ibm_marrakesh': (2, 4, 8), 'ibm_kingston': (2, 4)}
allr = []; far = []
for c in chips:
    for N in paperN[c]:
        cost = Rv[(c, 'R4', N)][0]; m, lo, hi = Rv[(c, 'R6', N)]
        rl, rh = abs(lo) / abs(cost), abs(hi) / abs(cost)
        allr += [rl, rh]; far.append(max(rl, rh))
        print('   %-14s N=%2d  no-echo cost R4 %+.5f | echo-minus-no-echo R6 %+.5f [%+.5f, %+.5f] | |end|/|cost|: lo %.2f, hi %.2f' % (c, N, cost, m, lo, hi, rl, rh))
print('   plain reading: each interval end divided by the no-echo cost (R4 median, [all]), N the paper says have no resolved change (2,4,8 fez/marrakesh; 2,4 kingston).')
print('   range over every end: %.2f to %.2f ; range of the larger end per case: %.2f to %.2f' % (min(allr), max(allr), min(far), max(far)))
print('   AMBIGUOUS: "interval ends reach 0.5 to 1.3 times the cost" has two plain readings; both are printed, neither is tuned.')
verdict(round(min(far), 1) == 0.5 and round(max(far), 1) == 1.3, '%.2f to %.2f (reading A: the farther end of each interval)' % (min(far), max(far)), '0.5 to 1.3')
verdict(round(min(allr), 1) == 0.5 and round(max(allr), 1) == 1.3, '%.2f to %.2f (reading B: every end, near ends included)' % (min(allr), max(allr)), '0.5 to 1.3')
