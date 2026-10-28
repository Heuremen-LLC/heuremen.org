# Addendum to PLAN_HAHN_FRESH.md and PLAN_HAHN_WIDER.md (22 Oct 2026 night), pushed BEFORE any fresh or wider run is read
After an outside referee re-read runs 1-2 from raw:
1. Runs of one chip are NOT combined by inverse variance (same 153 qubits, same stated values; per-qubit ratios correlate 0.989 kingston,
   0.853 marrakesh). They are pooled at the qubit level with hahn_pool.py. The combined values 0.922 [0.885, 0.959] (kingston) and
   0.933 [0.892, 0.975] (marrakesh) printed in results_hahn.text are WITHDRAWN and replaced by the pooled readings in pool_runs12.text.
2. Every chip is read with two fits: FIXED asymptote 1/2 (the committed hahn.py fit) and FREE asymptote (the model of IBM's open-source
   T2HahnAnalysis, base free). The fit-model choice moves the ratio by up to 7 % (fez 1.051 fixed, 0.979 free).
3. Bill rule, replacing the earlier one: a chip's line goes on the bill only if (a) its fresh-record run alone (stated T2 median age under
   3 h at firing) reads BELOW under BOTH fits, and (b) all its runs pooled by qubit read BELOW under BOTH fits. Anything else: site only,
   every run and both fits side by side.
4. The run-to-run shift under one unchanged IBM record (kingston 0.910 -> 0.946, marrakesh 0.907 -> 0.948, about 20 min apart) is
   larger than shot noise (referee: about 4.5 and 3.5 sigma) and is reported as delivered-T2 variation between jobs.
