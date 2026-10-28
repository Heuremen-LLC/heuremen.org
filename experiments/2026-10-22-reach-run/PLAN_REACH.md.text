# Reach run: where does the readout kick stop?
Written and SHA-256 hashed BEFORE any job is submitted (hash in PLAN_REACH.sha256; code hashes in CODE_REACH.sha256 and below).
Heuremen LLC (Wayfinder and Bones), 22 October 2026. Circuits: reach.py (dry run by default; nothing submitted). Analysis: analyze_reach.py,
which imports the deciding runs' ruler (deciding_runs/analyze.py) unchanged. Synthetic test: synth_reach.py (output synth_reach/SYNTH_RESULTS.txt).

## Question
On ibm_fez, ibm_kingston and ibm_marrakesh, each mid-circuit readout of 39 measured qubits shifts the phase of untouched qubits by about
1.6-2.9 deg per readout compared with an equal-length idle. The DISTANCE run (deciding_runs, run 1) read 0.81 / 1.77 / 1.30 deg at 3-4
couplings away. A referee split made after the data found the kick still present at exactly 5 couplings (8 qubits) and about 0 at 6 or more
(5 and 20 qubits). Where the kick stops is NOT established. This run reads the kick size at each exact distance d = 3, 4, 5, 6, 7 and 8 or more,
with about 38-42 untouched qubits in each bin on each chip.

## Design
- Distance d of an untouched qubit = number of couplings to the NEAREST measured qubit (breadth-first on the coupling map). Bins: 3, 4, 5, 6, 7, 8+.
- Three packed layouts per job. Each layout is a breadth-first ball of 39 measured qubits around a centre (ordered by distance, then qubit index).
  Rule for the centres (choose_layouts in reach.py; coupling map only, no data). Every qubit is a candidate centre. Every triple of centres is scored
  per bin on DISTINCT physical untouched qubits at that exact d in at least one layout. The search takes the triple that maximises, in order:
  (1) the smallest distinct count over d = 3-7, (2) the smallest distinct count over all six bins, (3) the smallest pooled (layout, qubit) count
  over d = 3-7, (4) the total pooled count. The last tie-break is a seeded random rank (seed 20261022).
  Balls centred in the middle of the chip win. A corner ball, such as the distance run's ball from qubit 0, has only 5-8 untouched qubits per exact
  shell. A central ball has 10-17.
- The three chips have the same 176-coupler map, so the rule picks the same centres on all three: 69, 71, 90.
- Per layout: idle (D) and readout (M) arms, N = 8 steps of length L (the chip's readout length, from qubit E[0]'s measure duration).
  In each step the 39 measured qubits are read (M) or idle for exactly L (D). The untouched qubits start on the equator (SX) and are never read
  mid-circuit. At the end they are read in X and in Y. The circuits are the control run's spectator_circuit, with a barrier fencing each step,
  identity layout, ALAP scheduling and optimization level 1. There is also one shared N = 0 baseline: all 156 qubits on the equator, read in X
  and Y. It feeds the coherence printout only.
- Pubs per copy: 2 baseline + 3 layouts x 2 arms x 2 settings = 14. Every circuit runs twice in its job (copy 0, copy 1), shuffled together
  (seed 20261029). Total 28 pubs, 4000 shots each.
- A unit is one untouched qubit in one layout (117 per layout, 351 per chip). The same physical qubit can sit in different bins in different
  layouts. The same physical qubit in the same bin in two layouts gives two units, so the distinct counts are printed beside.
- The distances used are the ones written into the job at build time. The analysis prints how many units would sit at a different distance on
  the couplers listed in that day's props.

## Layout counts from the dry run (identical on ibm_fez, ibm_kingston, ibm_marrakesh)
    untouched qubits at exact distance    d=1   d=2   d=3   d=4   d=5   d=6   d=7   d8+
    layout centre 69                       13    12    16    12    14    11    15    24
    layout centre 71                       11    10    14    11    15    10    12    34
    layout centre 90                        9    13    12    15    10    17    11    30
    pooled (layout, qubit) units           33    35    42    38    39    38    38    88
    distinct physical qubits               33    35    42    38    39    38    38    55
Readout length L: fez 415 dt, kingston 545 dt, marrakesh 671 dt (the same on every qubit of each chip).

## Seconds (dry run, 22 Oct 2026; nothing submitted)
ibm_fez 34 s (circuits 1.7-15.0 us); ibm_kingston 34 s (2.2-19.7 us); ibm_marrakesh 35 s (2.8-24.2 us). Total about 103 s. Remaining at the dry
run: 780 s.

## Expected widths (from the distance run's raw files, before this design)
In the distance run (4000 shots per copy), the per-qubit kick at d >= 6 scattered by a robust SD of 0.08-0.27 deg (median about 0.2). That is
about the shot noise. With about 38 units per bin, the median's 95 % half-width is about 1.96 x 1.25 x 0.2 / sqrt(38), roughly 0.08 deg per
chip, and about 0.05 deg combined. The d >= 8 bin (88 units) comes to about 0.05 per chip. That meets the 0.15 deg target wherever the kick is
small. Where the kick is present (d = 3-5), the distance run's per-qubit spread was REAL and large (robust SD 0.7-4 deg). More shots do not
narrow that. At d = 5, if the kick persists, expect a per-chip half-width of about 0.3-1.5 deg.

## The ruler (unchanged from the deciding runs)
Per unit: phase gap readout minus idle at N = 8, divided by 8, in deg (phase = atan2(<Y>, <X>), wrapped). Median over units, 95 % bootstrap
(2000 resamples). One of four outcomes: ABOVE (interval above 0), BELOW (below 0), CONSISTENT (contains 0; the interval is printed as the
bound), CANNOT TELL (fewer than 5 units, or the copies shifted).
Copies: per unit copy 0 - copy 1. If that interval contains 0, the pooled value is the reading. If it excludes 0, the reading is CANNOT TELL and
both copies are printed. Signs are carried: size = kick x the sign of the overall median (RR0). The sign of the d <= 2 median is printed beside
as a check. Chips are combined only by inverse variance (sigma = 95 % width / 3.92), with chi2/dof of the chips about the combination printed.
A chip whose copies shifted is listed and left out. Chips are never combined by counting labels.
DEAD-CLOCK GUARD: a unit whose coherence in the idle OR the readout arm at N = 8 (both copies pooled) is below 5 x shot noise
(5 / sqrt(8000) = 0.056) is left out of every phase reading. The count left out in each bin is printed. The N = 0 baseline median coherence is
printed beside.

## Readings and what each outcome means
- RR0 kick per readout, all units, all layouts (per chip; combined). It sets the sign. BELOW or ABOVE: a kick is present in this job (expected
  BELOW, as before). CONSISTENT: the overall median does not fix the sign. The sizes still carry it, and the d <= 2 sign check tells whether it
  agrees.
- RR1 kick size at d = 3, 4, 5, 6, 7, 8+ (per chip; combined).
  ABOVE: the kick is resolved at that exact distance. CONSISTENT: not resolved; the interval is the bound on its size there. BELOW: a kick of
  opposite sign to the overall one at that distance. CANNOT TELL: the copies shifted or there are fewer than 5 units.
- RR2 size(d) - size(d+1) for 3-4, 4-5, 5-6, 6-7, 7-8+ (per chip; combined; difference of medians, independent bootstrap of the two bins).
  ABOVE: the kick is smaller one coupling further out (it falls across that step). CONSISTENT: no resolved change across that step.
  BELOW: larger one step further out. CANNOT TELL: either bin's copies shifted, or the copy gap of the difference excludes 0.
- RR3 reach edge (per chip; combined). The largest bin b such that every bin from d = 3 up to b is ABOVE or CANNOT TELL, and b itself is
  ABOVE. A CANNOT TELL bin neither extends nor breaks the run. Printed beside: the first bin past the edge with its interval (the bound), and
  any ISOLATED ABOVE bins beyond the first CONSISTENT/BELOW bin, which are never taken as the edge.
  * Edge = 5 with 6, 7, 8+ CONSISTENT: the kick stops between 5 and 6 couplings, to within the printed bounds.
  * Edge = 6 or 7: it reaches that far. RR2 at the edge says whether it falls off sharply (ABOVE) or fades (CONSISTENT).
  * Edge = 8+: the kick reaches past 8 couplings. RR2 then reads the decay shape.
  * Edge = none: no bin resolved from d = 3. Either the kick is confined to d <= 2 in these layouts, or it cannot be told apart here.
  * Edge = 4 or 3 with d = 5 CONSISTENT: the referee split's d = 5 signal (8 qubits, one layout) does NOT repeat in 39 units.
  Before reading the combined edge, check chi2/dof: if it is far above 1, the chips disagree, and the per-chip edges are the reading.
- Context, not readings: units and distinct qubits per bin after the guard; per-layout median size per bin; the distance check against props.

## Synthetic test (synth_reach.py; FakeFez / FakeKingston / FakeMarrakesh; output synth_reach/SYNTH_RESULTS.txt)
The synthetic job, raw and props files use the real schema: raw bits are packed exactly as SamplerV2 returns them. The real layouts and
circuits list are used.
Realism planted:
- per-qubit coherence 0.45-0.80, with the readout arm 3 % lower;
- 3 dead clocks per chip;
- a static phase per qubit of -60..60 deg in both arms;
- a per-unit kick spread x exp(0.5 g).
The copies are identical in truth.
- (i) NO KICK: every RR0 / RR1 / RR2 reading on every chip and combined came out CONSISTENT; RR3 edge none on all chips and combined. The
  per-chip half-width at d = 6 was about 0.06 deg, and combined about 0.035 deg.
- (ii) STEP (kick -2.0 deg to d = 5, 0 beyond): RR3 edge = 5 on all three chips and combined. d = 6 combined read +0.019 [-0.016, +0.053]
  CONSISTENT; RR2 5-6 combined read +2.04 [+1.81, +2.26] ABOVE. All six combined bins covered the planted value. On fez, d = 3 came out
  CANNOT TELL (copy flag); the edge rule carried through it.
- (iii) SMOOTH (-3.0 exp(-(d-1)/2.5) deg): every bin ABOVE to d >= 8 on all chips and combined (edge 8+, "reaches past d = 8"). Every combined
  RR2 step came out ABOVE: falling at every step. Five of six combined bins covered the planted median. d >= 8 read +0.143 [+0.118, +0.168]
  against a planted +0.111. That is a small upward bias of the median in the widest bin, where the planted kick varies across d = 8-12 and the
  spread is skewed.

## FLOOR (six null seeds, no kick planted; same shots and analysis)
- Per chip: 13 of 216 readings resolved ABOVE/BELOW (6.0 %; 6 above, 7 below).
- Combined: 4 of 72 (5.6 %; 3 above, 1 below).
- RR3 named a reach edge on 0 of 18 chip-runs and 0 of 6 combined.
- An isolated ABOVE bin appeared beyond the (absent) edge in 2 of 6 combined null seeds. This is why the edge must be an unbroken run from d = 3.
- The copy flag turned 18 of 216 per-chip null readings into CANNOT TELL (8.3 %). RR2 inherits a flag from either bin, so it flags more often.
  The flag can only turn a reading into CANNOT TELL.

## Not established by this run, said beforehand
Distance on the coupling graph may not be the controlling variable. In the distance run, d = 5 sometimes read LARGER than d = 3-4 (kingston
3.1 vs 1.3-2.3 deg; marrakesh 2.3 vs 0.5-1.4). The per-qubit spread at d <= 5 was 0.7-4 deg. Something with a different geometry, such as shared
readout lines, could sit under a distance trend. This run reads distance bins only. The per-layout context lines show whether the bins agree
across the three centres, and they are not a reading.

## Code (sha256)
    27aebc5a2611dd73116c678a744a63c8bd1e64ec14373fd0aac224ec177b5a3c  reach.py
    65b6a0c56345d2b25de0d8c3002d11722c3be09006acb309685f10f5603e7cbb  analyze_reach.py
    362477139a0d3e75a2547dd03b82f90ad6feb3c03f5f65180cb4fbe0c025a7ad  synth_reach.py

Fire order (on his word only): python reach.py ibm_fez --fire; ibm_kingston; ibm_marrakesh. Then python analyze_reach.py > results_reach.txt.
