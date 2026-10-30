# Figures for "One Readout, Three Clocks"

Run from this folder with a Python that has numpy + matplotlib (built with `/opt/homebrew/bin/python3`, matplotlib 3.11.2, numpy 2.5.3;
`~/bones-forge/venv_quantum` has no matplotlib). Each script reads only files on disk and writes `fig<N>.pdf` (vector) and `fig<N>.png`.
`figstyle.py` holds the shared style (Okabe-Ito palette, 8-9 pt, 3.4 in / 7.0 in widths). Nothing in `main.tex` or the data folders is touched.

```
python3 fig1.py; python3 fig2.py; python3 fig3.py; python3 fig4.py
# fig5.py is a DRAFT: do not run into this folder until the reach result is written up (python3 fig5.py --out <dir> to preview)
```

**Sign convention.** Phase kicks are carried exactly as the analyses print them. The control-run, raw-readout, deciding-run and
kingston-repeat analyses all read the Y axis with the convention under which the kick is *negative*; the paper text quotes kick sizes
with the first run's positive sign (1.62, 2.86, 2.39 deg). Costs are idle minus readout (positive = less after readouts).
Open markers = the committed analysis printed "copies disagree" (CANNOT TELL); the pooled value and interval are still shown.

| Fig | Script | Width | Inputs (read only) |
|---|---|---|---|
| 1 | `fig1.py` | 3.4 in | `control_run/backend_properties_ibm_fez.json` (cz couplers), `~/bones-forge/heuremen.org/experiments/2026-10-04-three-positions/data/reads/countwide_roles_ibm_fez.json` (ROLES of `run.py`), `control_run/job_ibm_fez.json` (E, S, L_dt; asserted equal to roles) |
| 2 | `fig2.py` | 7.0 in | `control_run/job_<chip>.json`, `control_run/raw_<chip>.npz`, `control_run/results_prereg.txt` (R1 recomputed from raw with `analyze.py`'s bootstrap stream and asserted equal to the committed line) |
| 3 | `fig3.py` | 3.4 in | `control_run/results_prereg.txt` (R1, R2), `control_run/results_check_l.txt` (forward group), `deciding_runs/results_excitation.txt` (RJ.3), `excitation_kingston/results_exk.txt` (RJ.3) |
| 4 | `fig4.py` | 7.0 in | `control_run/results_iq.txt` (R1 clean, R2 all shots), `control_run/results_pairs.txt` (R1 clean, R3 all shots) |
| 5 (draft) | `fig5.py` | 3.4 in | `reach_run/results_reach.txt` (RR1 per distance bin) |

## Caption drafts

**Fig. 1.** (a) One step of each arm, on a shared time axis of exactly the readout length $L$ (415 dt = 1.660 µs on ibm_fez). Untouched
qubits are prepared on the equator ($\sqrt{X}$), then $N$ steps follow, each bounded by barriers on all qubits. In the idle arm the 39
measured qubits wait $L$; in the readout arm they are measured. Untouched qubits receive no gates during the steps in either arm and are
read out at the end in the $X$ or $Y$ basis. (b) Coupling map of ibm_fez (156 qubits, heavy-hex, cz couplers from the calibration at
firing) with the 39 measured qubits (squares) and the 117 untouched qubits (circles); the same roles are used on every processor.

**Fig. 2.** Phase of the 117 untouched qubits after $N$ steps, relative to $N=0$, in the idle arm and the readout arm of the same job
(control run; plan and analysis timestamped before the data). Points: median over qubits of the unwrapped per-qubit phase,
$\arctan(\langle Y\rangle/\langle X\rangle)$, both copies pooled; bands: 95% bootstrap intervals (2000 resamples). Text: the kick per
readout, (readout − idle phase at $N=8$)/8, median over qubits with 95% interval: fez −1.62° [−2.08, −1.18], kingston −2.86°
[−3.55, −2.09], marrakesh −2.39° [−3.00, −1.72] (negative in this analysis's $Y$ convention; sizes 1.62, 2.86, 2.39°).

**Fig. 3.** Controls on the kick. (a) Kick per readout with all 39 measured qubits read out, with 10 and with 1 read out on the same
schedule (the rest idling exactly $L$; control run, $N=8$), and for the minority of untouched qubits whose idle phase turns the other way
(21–36 per processor; control and raw-readout runs). Reading out 1 qubit leaves fez −0.08° [−0.15, +0.02], kingston +0.01° [−0.09, +0.08],
marrakesh −0.03° [−0.15, +0.10] (copies disagree); the minority is kicked the same way as the rest (fez −1.18, −0.96°; kingston −3.03,
−2.80°; marrakesh −3.37, −1.52°). (b) Kick with the measured qubits held in $|1\rangle$ minus kick with them in $|0\rangle$, each against an
idle in the same state: fez +0.18° [+0.02, +0.26], kingston +0.07° [−0.01, +0.23], marrakesh +0.92° [+0.53, +1.20] (deciding run), and the
kingston repeat +0.05° [−0.10, +0.13], shown side by side and not pooled. In this sign convention a positive value means a smaller kick in
$|1\rangle$. Medians over 117 untouched qubits, 95% bootstrap intervals; open markers: the two copies of the circuit disagree.

**Fig. 4.** Cost of the readouts to untouched qubits, idle minus readout, with all shots and with only shots in which every measured
neighbour was verified in its $|0\rangle$ cluster on every readout (97–99% of shots kept). (a–c) Coherence $|\langle X\rangle+i\langle Y\rangle|$
of the 86 untouched qubits wired to a measured qubit (raw-readout run), $N=2,4,8$; clean shots: fez +0.0042, +0.0057, +0.0043 (n.s.);
kingston +0.0054, +0.0056, +0.0062; marrakesh +0.0074, +0.0087, +0.0205. (d–f) Concurrence of untouched Bell pairs (wired, idle concurrence
above 0; pairs run), $N=2,4$; clean shots at $N=4$: fez +0.0070, kingston +0.0100, marrakesh +0.0093. Medians with 95% bootstrap
intervals, both copies pooled; open markers: the two copies disagree.

**Fig. 5 (draft, not final).** Size of the kick per readout against the number of couplings to the nearest measured qubit, with the 39
measured qubits packed in one region (reach run, three layouts per processor), median over untouched qubits per bin with 95% bootstrap
intervals. [Numbers to be filled once the reach write-up is done.]
