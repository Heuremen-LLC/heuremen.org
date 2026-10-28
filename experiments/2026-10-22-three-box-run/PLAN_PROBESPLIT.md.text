# Probe split: where does the three-box probe offset come from?
Heuremen LLC, Wayfinder and Bones, 22 October 2026 (night). Written, SHA-256 hashed (PLAN_PROBESPLIT.sha256) and pushed to heuremen.org
BEFORE any job of this run is submitted.

## What is known (not this run's evidence)
In the box whose ideal probe value is 0 (the unused box), the three-box probe reads +0.0812, +0.0617, +0.0485 (+- 0.010) on ibm_marrakesh
(22 Oct, variants V11, V10, VA; mean over all probe strengths, results_addendum_marrakesh.text R4); a standard noise model built from IBM's
calibration gives about +0.005, and the probe qubit's stated readout asymmetry about +0.006. The same offset appeared on ibm_fez (11 Oct,
+0.062). QIA-2026-001 row 12 carries it.

## Question (four answers allowed)
Is the offset (a) the probe's own X-readout bias, (b) what the system's state preparation and post-selection rotation do to the probe,
(c) what the controlled probe coupling adds, or several of these, or none resolved?

## Design (probe_split.py), same qubits as the three-box run
ibm_marrakesh 6, 7, 8 (s0, s1, probe; the 22 Oct layout); ibm_fez 142, 143, 144 (the 11 Oct and tonight's layout). Open instance of
.ibm_quantum_key_3. Five circuits, each twice (copy 0 / copy 1), shuffled (seed 20261023), 10000 shots, optimization level 1:
  K  probe |+> read in X;  P0  probe |0> read in X, system idle;  P2_V11, P2_V10, P2_VA  each variant's pre-selection and post-selection
  rotation with NO controlled gate, probe |0> read in X, kept shots s0 = s1 = 0.

## Readings (probe_split.py --analyze CHIP; bootstrap 4000; copies compared by their gap, a gap interval excluding 0 -> CANNOT_TELL)
  P0 <X>: the probe's readout bias.   P2 - P0: system operations.   S - P2: the coupling gate, where S is the unused-box offset of the
  same chip fixed before this run: marrakesh from results_addendum_marrakesh.text R4 (offsets_ibm_marrakesh.json, written before this
  plan); fez from tonight's three-box job dbd96jvi0auc7394vk40, R4 computed by the committed addendum_fourstate.py, unchanged.
  Each part is ABOVE / BELOW / CONSISTENT with 0 (the interval is the bound), or CANNOT_TELL when copies disagree.

KNOWN-GOOD: K must read <X> within 0.03 of 1 - 2 x IBM's stated prob_meas1_prep0 of the probe qubit (record saved at firing); if not met,
that chip's readings are not decided.

## Synthetic test (probe_split.py --synth WORLD -> synth_probe_split.text), run before this plan was hashed
Worlds null / readout (+0.06 in P0) / system (+0.06 in P2) / both, coupling part planted at +0.03 in all: readout and system parts are
recovered in every world (resolved when planted, consistent when not); the planted +0.03 coupling part resolves in 1 of 3 variants, so this
run names parts of about 0.05 and larger and leaves smaller ones as bounds.

## What this run cannot do
It does not test a coupling-gate contribution below about 0.05, does not separate readout bias from a bias of the final H gate, and the
marrakesh offset comes from a job 10 hours earlier on a different calibration cycle.
