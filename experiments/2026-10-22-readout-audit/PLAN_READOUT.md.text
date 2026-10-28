# Delivered readout error against IBM's stated readout fields, same hour (22 Oct 2026, night)
Written, hashed and pushed BEFORE any job. Code: readout_audit.py, fetch_ro.py (CODE_READOUT.sha256).
IBM's definition (quantum.cloud.ibm.com/docs/en/guides/qpu-information, read 22 Oct 2026): "The readout error quantifies the average
probability of incorrectly measuring a qubit's state. It is commonly calculated as the mean of prob_meas0_prep1 and prob_meas1_prep0";
"Readout error is measured by preparing the qubit in the |0> (|1>) state and measuring the probability of an output in the |1> (|0>)
state. The reported value is the average of these two errors."
Record finding that motivates the run (free, IBM's saved records, defined after reading them): in 154 of 5,387 distinct qubit measurements
(chip, qubit, measurement date) the stated readout_error differs from the mean of its own two components by more than 0.001, same date;
18 by more than 0.01, largest 0.157 (marrakesh). Stated component values are multiples of 1/4096; 397 qubit-records state P(1|0) = 0 exactly.
Design: ibm_fez, ibm_kingston, ibm_marrakesh (open instance of .ibm_quantum_key_3). P0: every qubit idle in |0>, measure all. P1: X on every
qubit, measure all. Each twice (copies), shuffled (seed 20261025), 10000 shots: 20000 shots per state per qubit. IBM record saved at firing.
Readings (readout_audit.py --analyze CHIP): R1 delivered (mean of P(1|0), P(0|1)) / stated readout_error, median over qubits, 4000-draw
bootstrap, four outcomes; same for each component. R2 on qubits whose stated readout_error is off its own definition by > 0.001: is the
delivered value nearer the stated readout_error or the stated mean? R3 qubits stated P(1|0) = 0: count with delivered P(1|0) Wilson 95 %
lower bound above 3/4096 (the 95 % upper bound of 0 of 4096). Copies compared by the median per-qubit gap.
KNOWN-GOOD: the median P(0) after |0> must be at least 0.90 per chip; if not, that chip is not decided.
Synthetic test (synth_ro.text): match 0.995 [0.988, 1.003] CONSISTENT; worse (1.5 x) 1.486 ABOVE; zeros_real: R3 finds 15 of 15, match finds 0 of 15.
Limits: all qubits are read at once (IBM's own simultaneity is not published); readout error drifts between IBM's cycles (record note 7),
so a ratio away from 1 is a same-hour reading, not a claim about IBM's average.
