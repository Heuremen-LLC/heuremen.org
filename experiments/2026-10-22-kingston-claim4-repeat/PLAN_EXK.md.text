# Claim 4, kingston repeat: is the readout kick bigger when the measured qubits are in 1?
Heuremen LLC, Wayfinder and Bones, 22 October 2026. Written and SHA-256 hashed (PLAN_EXK.sha256, code in CODE_EXK.sha256) BEFORE any job is submitted.

## Why
Claim 4 rests on the pre-registered excitation run (../deciding_runs/PLAN_EXCITATION.md). G = kick in 1 minus kick in 0, per step:
fez +0.178 [+0.025, +0.261], kingston +0.066 [-0.005, +0.229], marrakesh +0.923 [+0.533, +1.195]. An outside referee found the weak point:
on ibm_kingston the kick in 0 (RJ.1, M0 - D0) failed its copy check (copies shifted -0.099 [-0.178, -0.030] -> CANNOT TELL), so kingston's
G rests on a baseline whose two copies disagree. This run repeats the same job on ibm_kingston only, once, with more shots, and asks again.

## What is identical
- Circuits, roles and arms: ../deciding_runs/excitation.py build_excitation, imported and called unchanged (39 measured, 117 untouched spectators
  on the equator; D0, M0, D1, M1; N = 0, 4, 8; settings XX and YY; same roles file countwide_roles_ibm_kingston.json).
- Fire path: ../deciding_runs/runs.py main, executed from its own source unchanged except the two strings listed below: every circuit twice
  (copy 0, copy 1), shuffled with seed 20261022 + 7 (checked: the same order as the original kingston job, 40 pubs), transpile optimization 1,
  identity layout, alap, seed_transpiler 7, Darkwave Disco instance, SamplerV2, the same ESTIMATED formula and the same seconds stop.
- Analysis: ../deciding_runs/analyze_excitation.py and analyze.py, unchanged, run as
  `python ../deciding_runs/analyze_excitation.py /Users/heuremen/bones-forge/arxiv_one_readout/excitation_kingston`
  (checked on a synthetic kingston folder: it reads the new files as they are). Same four readings, same rulers (median over spectators,
  2000-resample bootstrap 95 %, copy check first, four outcomes).
- Those files' hashes equal ../deciding_runs/CODE_EXCITATION.sha256; they are listed again in CODE_EXK.sha256.

## What changed (all in excitation_k.py, nothing else)
1. Shots per circuit 2000 -> 4900 (the largest round count with the dry-run ESTIMATED at or under 60 s).
2. Chip: ibm_kingston only.
3. Output folder: this folder (runs.py's HERE set to it). File names unchanged (job/props/raw_excitation_ibm_kingston.*), so the analysis is unchanged.
4. Plan hash written into the job file: PLAN_EXK.sha256 instead of PLAN.sha256 (runs.py) / PLAN_EXCITATION.sha256 (excitation.py).
5. Two extra keys in the job file's info record: repeat_of = dbcnds424mrc73ck0450 (the original kingston job), shots_original = 2000.
   The analysis reads only info['S'].
Not changed on purpose: the shuffle seed. It puts copy 0 slightly earlier in the job than copy 1 (mean position 17.4 vs 21.6 of 0-39);
same as the original, so the copy check means the same thing in both jobs.

## Readings (unchanged from PLAN_EXCITATION.md), kingston only
RJ.1 kick in 0 (M0 - D0); RJ.2 kick in 1 (M1 - D1); RJ.3 G = RJ.2 - RJ.1; RJ.4 measured qubits in 1 without readout (D1 - D0). deg/step.
Sign: the kick is read with its sign. If RJ.1 and RJ.2 are both BELOW (kick negative, as on all three chips on 21 Oct), RJ.3 ABOVE means the
kick is SMALLER in 1 and RJ.3 BELOW means BIGGER in 1; if both are ABOVE, the reverse; if their signs differ, RJ.3 is reported as values only.
(PLAN_EXCITATION.md's line "ABOVE -> bigger in 1" assumed a positive kick; this paragraph states it by size.)

## What each outcome means for claim 4 on kingston
- RJ.1 and RJ.3 copies agree, RJ.3 says SMALLER in 1: kingston joins fez and marrakesh on a baseline that passed its copy check.
- RJ.1 and RJ.3 copies agree, RJ.3 CONSISTENT: no difference resolved on kingston; claim 4 on kingston reads "not bigger by more than the
  interval's bigger-in-1 end", and that end is the number quoted.
- RJ.1 and RJ.3 copies agree, RJ.3 says BIGGER in 1: claim 4 does not hold on kingston; it is narrowed to chip-dependent, fez and marrakesh.
- RJ.1 copies shift again (CANNOT TELL): kingston's kick-in-0 baseline disagrees with itself in two jobs out of two; kingston is reported as
  CANNOT TELL for claim 4 and the claim cites fez and marrakesh only. More shots make this check stricter, not looser: the copy-gap interval
  is about 0.65x as wide as on 21 Oct, so a real shift of 0.1 deg/step would show; a chance flag (the 21 Oct interval ended at -0.030, just off zero)
  passes about 95 % of the time.
- RJ.3 copies shift: CANNOT TELL on kingston, whatever RJ.1 does.
- The 21 Oct kingston values are printed beside the new ones; the two jobs are NOT combined (same 117 spectators, not independent qubits).
A resolved |G| under about 0.04 deg/step with an interval just off zero sits at the null floor (below).

## Synthetic test (synth_exk.py: synth_excitation.py's own source on FakeKingston, 117 spectators, contrast 0.8, 6 seeds)
4900 shots (results_synth_exk_4900.txt):
- null G = 0: RJ.3 read +0.004 mean (-0.030 to +0.041), interval width 0.071 mean (0.050 to 0.086); 5 CONSISTENT, 1 BELOW (seed 1,
  -0.030 [-0.061, -0.002]: the 1-in-20 a 95 % interval gives); copy flags on RJ.3 0 of 6; across the four readings 1 flag in 24.
- planted G = +0.5: RJ.3 read +0.497 mean (+0.466 to +0.533), width 0.067; 6 of 6 ABOVE, every interval contains +0.5.
2000 shots (results_synth_exk_2000.txt, same seeds): null RJ.3 width 0.110 mean; planted read +0.485, width 0.117; null 6 of 6 CONSISTENT.
The synthetic gives every spectator the same G, so 0.071 is the SHOT-NOISE floor of the width, not the expected real width.

## Expected interval width on the chip
The 21 Oct kingston RJ.3 interval is 0.234 wide. Its per-spectator spread splits (copy-gap method, robust sd) into shot noise 0.31 and
spectator-to-spectator spread 0.52 deg/step: about a quarter of the variance is shots. At 4900 shots the shot part falls from about 0.110 to
0.071 in width; the spread does not move. Expected RJ.3 width about 0.22 (1.07 to 1.09x narrower), and no shot count goes below about 0.21.
A 2x narrower interval is NOT reachable by shots on one chip with this ruler; it would need about 4x as many independent spectators.
What this run buys is an independent second kingston job with a fresh copy check on RJ.1, at a slightly narrower interval.

## Dry run (dryrun_exk_ibm_kingston.log, 22 Oct 2026)
excitation ibm_kingston: 40 pubs, 4900 shots each, circuits 2.2-19.7 us, total shots 196000, ESTIMATED 59 s (account had 780 s left, queue 1).
Fire: `cd /Users/heuremen/bones-forge/arxiv_one_readout/excitation_kingston && ~/bones-forge/venv_quantum/bin/python excitation_k.py --fire > fire_exk_ibm_kingston.log 2>&1`
