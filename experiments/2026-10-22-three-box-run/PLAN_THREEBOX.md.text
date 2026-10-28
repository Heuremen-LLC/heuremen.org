# Three-box weak values: does the shortfall follow the negative sign or the corner?
Heuremen LLC, Wayfinder and Bones, 22 October 2026. Written and SHA-256 hashed (PLAN_THREEBOX.sha256; code in CODE_THREEBOX.sha256) and pushed
to heuremen.org BEFORE any job of this run is submitted.

## What is already known (read before this plan; not part of this run's evidence)
- 11 Oct 2026, ibm_fez, qubits 142/143/144, 9 jobs, probe strength eps 0.1-0.4 (circuits read back from IBM; generating script not in our
  record): fit <Xp> = S + P sin(2 eps) gives the negative ("corner", |11>) box P = -0.530 +- 0.053 (textbook -1), positive boxes +0.817 and
  +0.819, unused box P -0.101 +- 0.060 with S +0.111 +- 0.031.
- 24 Sep 2026, eps 0.25: negative box fez -0.498, kingston -0.587, marrakesh -0.689 (standard variant, negative at |10>), marrakesh -0.506
  (corner variant). IBM's noise model (FakeMarrakesh) gave -0.81 for the negative box.

## Question
The negative box reads about half of textbook; the positive boxes about 0.8. Does the shortfall follow the NEGATIVE SIGN (whichever box
carries it) or the POSITION (the |11> corner)? And how does each box compare with IBM's own noise model built from the calibration in force
when the job fires?

## Design (tb_common.py, run_threebox.py)
- Two system qubits (s0, s1) and one probe. Pre-select psi = uniform over three basis states; weakly rotate the probe RY(2 eps) only if
  the system is in the probed box; rotate the post-selected state phi to |00>; read the probe in X. Kept shots: s0 = s1 = 0.
- Three variants, each probing all four boxes:
  V11: support {00, 01, 11}, phi = (00 + 01 - 11)/sqrt3   (negative box at the corner, as 24 Sep / 11 Oct)
  V10: support {00, 01, 10}, phi = (00 + 01 - 10)/sqrt3   (standard Aharonov-Vaidman; negative box at 10; corner unused)
  VA : support {00, 01, 11}, phi = (-00 + 01 + 11)/sqrt3  (same support as V11; negative box moved to 00; corner now POSITIVE)
- eps = 0.10, 0.15, 0.20, 0.25, 0.30, 0.40 (the 11 Oct grid). 3 variants x 4 boxes x 6 eps = 72 circuits, each twice (two copies),
  shuffled in one job per chip (seed 20261022). State preparation by qiskit StatePreparation, transpiled at optimization level 1.
- Chips and budget (seconds remaining at writing: open instance on .ibm_quantum_key 357 s; Darkwave Disco 271 s):
  ibm_fez, qubits s0,s1,probe = 142, 143, 144 (the 11 Oct qubits), 7000 shots per circuit, open instance (.ibm_quantum_key);
  ibm_marrakesh, layout by the fixed rule in run_threebox.choose_layout (lowest readout + two-qubit error chain at firing; dry run gave
  6, 7, 8), 6000 shots per circuit, Darkwave Disco. Kingston is not run (no budget).
- IBM's prediction: at firing the backend properties are saved, and IBM's noise model (qiskit_aer NoiseModel.from_backend on that backend)
  is evaluated on every transpiled circuit by density-matrix simulation and saved (model_<chip>.json) BEFORE the job result is read.
  The noiseless value of every logical circuit is saved beside it.

## Readings (analyze_threebox.py; every one per chip; four outcomes; copies compared by their difference, a gap excluding 0 = CANNOT_TELL)
- For every variant and box: P from the fit <Xp> = S + P sin(2 eps) over the six eps (copies pooled), 95% interval by binomial
  resampling within each circuit (2000), with the same fit applied to the IBM-model and noiseless values.
- R1 (decision): negative box, hardware P minus IBM-model P, per variant.
- R2 (decision, the question): G = (hardware P / IBM-model P on the negative box) minus (mean of that ratio on the two positive boxes),
  per variant. If the shortfall follows the SIGN, G is BELOW 0 in all three variants. If it follows the CORNER, G is BELOW in V11,
  about 0 in V10 and ABOVE in VA (its corner box is positive). Any other pattern is reported as it falls.
- R3: the unused box's S offset and P, per variant.
- Across chips: inverse-variance combination of G per variant, chi2/dof printed, combined value quoted only if chi2/dof < 2.
- Everything else is defined after the data and labelled so.

KNOWN-GOOD: the committed fit and kept-shot reader (known_good.py) applied to the 11 Oct ibm_fez three-box jobs must reproduce corner P = -0.530 and box-A P = +0.817 to within 0.01 -- it does (known_good_synth_check.text); and on hardware in this run the noiseless-value fit of the positive boxes must lie within 0.1 of +1 (finite-eps check).

## Synthetic test (synth_threebox.py -> synth_results.text), run before this plan was hashed
Three planted worlds on three synthetic chips at 6000 shots: NULL (hardware = model) reads CONSISTENT (one variant not combined, chi2/dof
2.11); SIGN (negative box x0.55) reads G BELOW in all three variants; CORNER (|11> box x0.55) reads BELOW / CONSISTENT / ABOVE for
V11 / V10 / VA. This run has two chips, not three, so its combined intervals are wider than the synthetic ones.

## What this run cannot do
- Two chips, one day, one layout per chip. eps = 0.4 is not in the weak limit; P is compared with the model and noiseless fits over the
  same grid, like for like. The IBM noise model is IBM's published calibration through qiskit-aer's standard construction; a mismatch
  measures that model, not the chip's full error budget.
