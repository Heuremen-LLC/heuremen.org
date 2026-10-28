# Direct per-edge zz: does IBM's stated zz predict the frequency shift an untouched qubit feels when one neighbour goes from 0 to 1?
Heuremen LLC, Wayfinder and Bones, 22 October 2026. Written, SHA-256 hashed (PLAN_DIRECT.sha256; code CODE_DIRECT.sha256) and pushed to
heuremen.org BEFORE any job of this run is submitted.

## Why
The pre-registered flat-beat run on ibm_phoenix (zz_flat, 37 fresh pairs) found the beat of an untouched qubit beside one SUPERPOSED
neighbour scales with IBM's stated zz with log-log slope 0.19 [0.01, 0.38] (ibm_fez, same fit, 1.00). An outside referee could not break it
(a simulation of the exact grid returns slope 1.00-1.02 for beats at the stated values) and named the cheapest test that removes the beat
fit entirely: Ramsey on A with the neighbour in 0 against in 1 -- the frequency shift is the zz itself, in IBM's own sense, with no envelope
model and no superposed neighbour. An attempt to read it from the 22 Oct one-neighbour run's |1> arm was judged not evidence (5 us steps,
drift), so this run measures it directly.

## Design (direct.py)
Triples (B1, A, B2) fixed before the job: ibm_phoenix = the 20 triples of the flat-beat job; ibm_fez = the 33 triples of the 22 Oct
one-neighbour job. A: SX, wait tau, read X or Y. Arms Z (B1, B2 in 0), O1 (B1 in 1), O2 (B2 in 1). tau = 1-39 us in 2 us steps (at the
largest stated phoenix zz, 59 kHz, 0.74 rad per step: no aliasing). Every circuit twice (copy 0 / copy 1), shuffled; 1000 shots.
ibm_fez on the open instance of .ibm_quantum_key; ibm_phoenix on Darkwave Disco when its usage cap allows (phoenix is not on the open
instance). IBM's calibration saved at firing (props_direct_<chip>.json) gives the stated zz.

## Readings (analyze_direct.py; unit = edge; bootstrap over triples, 2000; copy gap excluding 0 -> CANNOT_TELL)
Per edge, heralded on the neighbour's final readout (O shots with it read 1, Z shots with it read 0), df = slope of
unwrap(phase_O - phase_Z) against tau / 2 pi.
D2 (DECISION) log-log slope of |df| on stated |zz| (edges stated > 1 kHz). TRACKS if the interval contains 0.99 (the synthetic tracking
value); DOES NOT TRACK if the interval lies below it; a factor-2 convention still reads ~1.0 here.
D3 median |df| / stated (a constant convention factor shows here, not in D2). D4 rank correlation. D5 edges IBM states as 0.00.
KNOWN-GOOD: on ibm_fez D2's interval must contain 0.99 and D3 must lie in [0.4, 2.5]; otherwise phoenix's reading is not decided.

## Synthetic test (synth_direct.py -> synth_direct_results.text), run before this plan was hashed
Phoenix triples, stated zz, A detuning, T2 60 us, neighbour T1 120 us, 1000 shots: TRACK reads D2 0.991 [0.985, 0.999], D3 0.999;
HALF (factor 2) reads D2 0.996, D3 0.500; FLAT 13 kHz reads D2 0.007, D3 0.459. The first version of the analysis read 5% of every planted
shift (1 - 2a computed on unsigned integers); the synthetic test caught it and it was fixed before this plan was hashed.

## What this run cannot do
One job per chip, one day; the neighbour's other neighbours sit in 0; a shift that depends on more than the one neighbour is not tested.
