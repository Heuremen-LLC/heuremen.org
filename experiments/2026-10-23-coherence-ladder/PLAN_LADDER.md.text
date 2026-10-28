# The coherence ladder: T1, T2*, Hahn T2, CPMG-4, CPMG-16 on the same qubits in one job (kingston, marrakesh)
Pre-registered. Written and SHA-256 hashed (PLAN_LADDER.sha256; code in CODE_LADDER.sha256) BEFORE any job is submitted. Heuremen LLC,
Wayfinder and Bones, 23 October 2026. Circuits, fire path and analysis: ladder.py (one file). Synthetic test: synth_ladder.text.

## Why this run
The Hahn T2 runs (../hahn_t2/) found IBM's published per-qubit T2 differs from the T2 we measure by about 20 to 60 % per qubit, with no
detectable shrink for a record under an hour old, no repeat per qubit across records, and no less than the drift between our own
measurements hours apart (precision.py, refereed 23 Oct). IBM publishes one echo T2 and one T1 per qubit. This run reads the same qubits
on five rungs in one job, so each rung is a second and third seat on the same qubit at the same moment, and the rungs together say what
kind of noise the qubit is drifting with: decoupling with more pulses filters out slower noise.

## Design
ibm_kingston and ibm_marrakesh, Darkwave Disco instance of .ibm_quantum_key_3. One sparse set (set 0 of hahn.py's sparse_sets: 92 qubits,
no two coupled); every other qubit idle in |0>. Per rung, 8 waits (us):
  T1      X, wait, measure                                        [0, 20, 50, 100, 170, 260, 400, 600]
  RAMSEY  SX, wait, (RZ(pi/2) or not), SXdg, measure               [0, 5, 10, 20, 35, 60, 100, 160]  two quadratures; coherence is
          the magnitude of the two (shot-noise floor removed), so a frequency offset cannot fake a decay
  HAHN    SX, wait/2, X, wait/2, SXdg (hahn.py's circuit)          [0, 20, 40, 90, 160, 220, 300, 400]
  CPMG-4  SX, then 4 Y pulses spaced wait/4 (half spaces at ends), SXdg    [0, 40, 90, 160, 260, 400, 600, 800]
  CPMG-16 the same with 16 Y pulses                                [same]
Barriers around every pulse; ladder.py refuses to fire if any pi pulse is missing after transpilation (dry runs: 0 failures; one CPMG-16
circuit carries 1,472 X pulses = 16 x 92). Every circuit twice (copy 0 / copy 1), all shuffled together, 1,200 shots. 96 circuits per chip,
estimated 54 s each (dry runs dryrun_ibm_*.log). IBM's record saved at submission (props_ladder_<chip>.json).

## Fits
T1: P1 = B + A exp(-t/T1), B free. RAMSEY: coherence c = A exp(-t/T2*). Echo rungs: both fits of the Hahn work, FIXED (asymptote 1/2) and
FREE (asymptote free); every reading printed under both. A qubit's rung counts if A > 0.4 (T1) or 0.2 (others) and T < 1990 us.

## Readings (median over qubits of a per-qubit quantity, 4000-draw bootstrap over qubits, four outcomes)
  R1  delivered T1 / IBM stated T1 (same record)                    vs 1
  R2  delivered Hahn T2 / IBM stated T2 (same record)               vs 1
  R3  T2(rung) / T2(HAHN) per qubit, rung = RAMSEY, CPMG4, CPMG16   vs 1
  R4  T2(rung) / (2 T1) per qubit, every echo rung                  vs 1 (the physical limit)
  COPY RULE: per rung, median T(copy 0)/T(copy 1); a rung whose interval clears 1 by more than 5 % is SHIFTED and every reading using
  it is CANNOT_TELL. (Synthetic floor with no planted shift: 3.2 %; a planted 20 % shift is caught on every echo rung.)
  NOISE CLASS (from R3 CPMG16 and R4), one of: T1-LIMITED (R4 Hahn upper bound >= 0.9) / SLOW NOISE (CPMG16 ABOVE Hahn, still BELOW
  2 T1) / SLOW NOISE, REMOVED (CPMG16 ABOVE, reaches 2 T1) / WHITE (CPMG16 CONSISTENT with Hahn) / decoupling SHORTENS T2 / CANNOT_TELL.

KNOWN-GOOD: two real-hardware readings whose answer is known, read in this job on each chip. (1) median Hahn P(1) at wait 0 at least 0.90
and median Ramsey coherence at wait 0 at least 0.85 (the hahn.py known-good, met on every earlier Hahn run). (2) Physics: no T2 can exceed
2 T1, so the upper 95 % bound of median T2(HAHN)/(2 T1) must be at most 1.05 under both fits. A chip that fails either has its readings
printed and not used.

## Bill rule (written now)
R1 or R2 BELOW under BOTH fits on a chip with known-good met: that line joins the T2 rule already in force (../hahn_t2/PLAN_HAHN_FRESH.md):
it goes on the bill only after a run against a fresh record reads BELOW alone and all runs combine BELOW (chi2/dof < 2). This run alone
puts no line on the bill. R3, R4 and the noise class are science (heuremen.org), not bill lines.

## Synthetic test (synth_ladder.text, 1,200 shots, 90 qubits, both fits)
Planted worlds and what they must read: slow (1/f-like) -> SLOW NOISE; white -> WHITE; t1limit -> T1-LIMITED; mixed (half slow, half
white) -> SLOW NOISE; short (every T1 and T2 at 0.85 x stated) -> R1 and R2 BELOW; pulse (0.4 % contrast lost per pulse) -> WHITE, not
read as a T2 change; drift (copy 1 echo T2 x 0.8) -> CANNOT_TELL. All seven read as planted under both fits.
