# Readout: repeat and isolated read of IBM's stated zeros (22 Oct 2026 night), pushed BEFORE these jobs
Run 1 (results_ro_*.text, zeros_history.text): 17 qubits stated P(1|0) = 0 (0 of 4096); same evening, all-qubit read, 20000 shots: 16 above the
0-of-4096 bound (marrakesh 14 of 14, median 0.0073; zeros measured 3.9 h before our job on marrakesh, 3.2 h fez, 0.9 h kingston).
The vendor's remaining answers: (a) the qubits drifted between IBM's measurement and ours; (b) reading every qubit at once raises P(1|0).
(a) REPEAT: readout_audit.py ibm_marrakesh --fire, unchanged. If IBM's record has not changed and delivered P(1|0) on the zero qubits is
    unchanged within its interval, the delivered side is stable over the interval; IBM's zero stays as stated. Run 1 files kept as *_run1.
(b) ISOLATED: readout_isolated.py CHIP --fire on fez, kingston, marrakesh: only the qubits stated zero AT FIRING are measured; everything
    else untouched. Reading: per-qubit isolated - all-qubit (run 1) P(1|0), median, bootstrap, four outcomes; and R3 on the isolated read.
Bill rule for the zero line: it stands as a deviation if the isolated read alone puts at least half of the stated-zero qubits above the
0-of-4096 bound, on the record the zeros were stated in. Otherwise it is a record note with every number side by side.
KNOWN-GOOD: readout_audit known-good (median P(0 | |0>) >= 0.90) for the repeat; for the isolated read, mean P(1 | |1>) over the listed
qubits must be at least 0.80. Synthetic test: synth_iso.text (same -> CONSISTENT; crosstalk -> BELOW) and synth_ro.text (unchanged).
