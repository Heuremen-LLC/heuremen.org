# Known-good seat for the flat-beat analysis, on real hardware (written and pushed BEFORE the flat-beat job's output is opened)
Heuremen LLC (Wayfinder and Bones), 22 October 2026. Added after the plan was pushed and the job fired, and before its data returned.
Why: analyze_flat.py was tested only on synthetic data. If phoenix reads FLAT, nothing in the run separates "phoenix's beat does not follow
the calibration" from "this free fit, on real hardware, flattens every beat". A seat where the answer is already known is needed.
The seat: knowngood.py runs analyze_flat's own fit and readings, unchanged, on the 22 Oct one-neighbour run (ibm_fez: an outside referee's
separate free fit found the beat tracks the calibration, slope 1.08 [0.94, 1.25]; ibm_phoenix: the same referee found it flat).
Rule, fixed now:
- If on ibm_fez the slope is ABOVE 0 and its interval reaches 1, the fit can see a tracking beat on real hardware: a FLAT reading on the new
  phoenix pairs is read as a property of phoenix.
- If on ibm_fez the slope is not ABOVE 0, the fit cannot see a tracking beat on real hardware on that grid, and a FLAT reading on the new
  phoenix pairs is CANNOT TELL (instrument), not a finding.
- The ibm_phoenix line is printed beside it and is not a test (it is the data the hypothesis came from).
Caveat stated now: the 22 Oct fez grid ran to 303 us and the new phoenix grid to 162 us, so this seat checks the fit, not the new grid.
