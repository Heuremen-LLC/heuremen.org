# Addendum to PLAN_THREEBOX.md -- written after the jobs were submitted and BEFORE any result was read (22 Oct 2026)

1. The plan's question "does the shortfall follow the negative sign OR the corner" has four answers, not two: SIGN, CORNER, BOTH, NEITHER
   (plus cannot tell). addendum_fourstate.py reads all of them per chip: r_b = hardware P / model P over the nine support boxes, regressed on
   [negative sign], [|11> corner], [|00> position]; the state is named from the sign and corner terms; copy 0 and copy 1 are compared on those
   terms. Planted test (synth_addendum_results.text, two synthetic chips per world, 6000 shots): SIGN, CORNER, BOTH and |00>-position worlds
   are named correctly on both chips; the NULL world reads NEITHER on one chip and SIGN on the other (a false SIGN at this shot count; so a
   SIGN on one chip alone is not a finding, and the two chips are read side by side).
2. Wording: "IBM's noise model" in the plan means qiskit-aer's standard NoiseModel.from_backend built from IBM's calibration at firing. It is
   a construction from IBM's stated figures, not a figure IBM states.
3. An outside referee of the 11 Oct fez data (REFEREE_11OCT_before_results.md, read before any result of this run) found: the raw numbers
   reproduce; a FakeFez noise simulation plus ONE common probe offset (+0.05) and gain (0.96) reproduces all 24 (box, eps) means (chi2 41.6 for
   22); the corner/positive asymmetry is already present in the noise model; what survives is a positive probe offset of +0.062 +- 0.010 in
   every box (24 Sep: +0.060 +- 0.016) that the noise model lacks. So this run's readings carry one more, defined now:
   R4: the common probe offset, mean <Xp> on kept shots of the unused box over all eps (hardware against model), per variant and chip.
   And R1/R2 are read with that in mind: a shortfall that disappears once a common offset and gain are allowed is not box-specific.
