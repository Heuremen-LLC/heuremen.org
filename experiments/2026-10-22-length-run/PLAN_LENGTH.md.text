# Two readout lengths on one chip: does the kick come per readout or per microsecond of readout?
Heuremen LLC, Wayfinder and Bones, 22 October 2026. Written, SHA-256 hashed (PLAN_LENGTH.sha256; code in CODE_LENGTH.sha256) and pushed to
heuremen.org BEFORE any job of this run is submitted.

## Why
The paper leaves two things open that this run can read: (1) the mechanism of the phase kick on untouched qubits -- a control-system
phase-frame update applied once per readout, or a physical frequency shift during the readout pulse, which would grow with its length; the
21 October run had only three readout lengths, one per processor, so a per-readout and a per-length kick could not be told apart; and
(2) the measured qubit's loss at equal elapsed time exists for ibm_fez only.

## Design (run_len.py)
- Same roles as the control run (39 measured, 117 untouched; countwide_roles_<chip>.json). Every readout step is the processor's `measure`
  (L1) or its second readout instruction `measure_2` (L2), each against an idle of exactly its own length; ibm_kingston 2.180 / 1.380 us,
  ibm_marrakesh 2.684 / 1.884 us, ibm_fez 1.660 / 1.340 us (from the backend target at the dry run).
- Arm S: the 117 untouched qubits on the equator (sx), N = 4 and 8 readouts or idles of the 39, read X and Y.
- Arm B: the 39 measured qubits prepared in 1, N = 4 and 8 readouts (L1 or L2) or idles of equal length, then every measured qubit read
  with `measure`.
- 24 circuits, each twice in the job (copy 0 / copy 1), shuffled (seed 20261023); 5000 shots per circuit; one job per processor on the
  open instance of .ibm_quantum_key. Transpiled at optimization level 1 with the identity layout and ALAP scheduling.

## Readings (analyze_len.py; four outcomes; a copy-gap interval excluding 0 makes a reading CANNOT_TELL)
- R1: kick per readout at L1 and at L2, N = 4 and 8: wrap(phase_readout - phase_idle)/N in degrees per untouched qubit, median over
  qubits alive in both arms (coherence > 5/sqrt(shots)), 95% bootstrap over qubits (2000).
- R2 (decision): ratio k(L2)/k(L1), same qubits, bootstrap. PER-READOUT if its interval contains 1 and not L2/L1; PROPORTIONAL-TO-LENGTH
  if it contains L2/L1 and not 1; CANNOT_TELL if it contains both (fez, where L2/L1 = 0.81, may land here); NEITHER (between or outside)
  otherwise, read through R3.
- R3: the split k = a + b L over the two lengths: a = per-readout part (deg), b = per-microsecond part (deg/us), with bootstrap intervals.
- R4: the measured qubit, P(1) after N idles minus after N readouts of equal length, median over the 39, at L1 and L2.
- Processors are reported side by side; the ratio is not combined across processors.

KNOWN-GOOD: on each processor the L1 (`measure`) kick at N = 8 must read BELOW 0 with an interval overlapping the control run's committed kick on that processor (same roles, same instruction, same N; control run: fez -1.62 [-2.08, -1.18], kingston -2.86 [-3.55, -2.09], marrakesh -2.39 [-3.00, -1.72] deg per readout). If it does not, R2-R4 on that processor are reported as not decided.

## Synthetic test (synth_len.py -> synth_len_results.text), run before this plan was hashed
Planted worlds on kingston-like and marrakesh-like lengths at 5000 shots: PER-READOUT reads PER-READOUT (ratio 0.97-1.02); PROPORTIONAL
reads PROPORTIONAL-TO-LENGTH (0.61-0.75); MIXED (a = -1.0 deg, b = -0.6 deg/us) reads NEITHER: between, and the split recovers
a -1.07 to -1.22, b -0.54 to -0.58; NULL reads CANNOT_TELL. One R4 copy gap fired falsely on a planted world (the floor of that rule).

## What this run cannot do
One layout and one day per processor; the measured qubits' roles are those of the control run. A per-readout kick is what a phase-frame
update would give, but any effect that happens once per readout regardless of its length (for example at the pulse edges) would read the
same; a per-length kick is what a frequency shift during the pulse would give.
