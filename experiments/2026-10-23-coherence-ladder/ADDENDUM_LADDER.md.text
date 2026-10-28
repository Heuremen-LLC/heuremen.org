# Addendum to PLAN_LADDER.md -- written after both jobs were fetched, BEFORE any result was read (23 Oct 2026, /q pass)
Both jobs ran within 2 s of submission (run_ladder_<chip>.json: kingston created 20:37:47Z, running 20:37:49Z; marrakesh 20:38:01Z /
20:38:02Z), so the record saved at submission is the record they ran under. IBM charged 52 s and 53 s. ladder.py is unchanged; every
registered reading prints from it as written. ladder_addendum.py adds the readings below (its docstring is the full rule):
  X1  T2(CPMG16)/T2(CPMG4): both Y-pulse rungs, so only the pulse count changes (R3 against HAHN also changes the pulse axis X -> Y).
      The NOISE class reads SLOW only if X1 is not BELOW; otherwise CANNOT_TELL (axis-confounded).
  X2  Ramsey refit on waits with coherence >= 5x its shot-noise floor; if its outcome differs from the registered R3 RAMSEY, RAMSEY is
      CANNOT_TELL.
  X3  dropped and capped qubits per rung; CPMG known-good: median 1-P(1) at wait 0 >= 0.90, a rung failing it is CANNOT_TELL.
  X4  T1 shape, descriptive: 1/e crossing over fitted T1. On a pure single exponential the interpolation itself reads 1.014 to 1.021
      (synth_ladder_addendum.text), so that, not 1, is the reference.
  X5  R1 and R2 against the record as it stood after the job, and the count of qubits whose T1/T2 changed in between.
  X6  transpiled duration minus nominal wait per rung (sched_ladder.text): constant inside each rung to 0.03 us; pulse and readout
      overhead 2.2 to 3.3 us, the same at every wait, so it shifts no fitted time.
  BOTH FITS: a NOISE class is reported only if FIXED and FREE give the same class; otherwise CANNOT_TELL. (Found in the synthetic test:
      a planted Gaussian-decay world reads SLOW under FIXED and T1-LIMITED under FREE.)
  TWO CHIPS: an R3, X1 or NOISE outcome is a finding only if both chips read it; one chip is printed as one chip.
Synthetic floor (synth_ladder_addendum.text): white and slow worlds, seeds 11 to 15, both fits: every class read as planted.
