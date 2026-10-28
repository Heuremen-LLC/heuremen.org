# IBM's QAOA tutorial, run as written, against the forecast from IBM's published figures
Pre-registered; hashed with the code (PLAN_TUTORIAL_QAOA.sha256, CODE_TUTORIAL_QAOA.sha256) BEFORE any job. Heuremen LLC, 23 Oct 2026.
Source procedure: quantum.cloud.ibm.com/docs/en/tutorials/quantum-approximate-optimization-algorithm (small example; settings quoted in
tutorial_qaoa.py's docstring). Chips: ibm_kingston and ibm_pittsburgh (the chip the tutorial's own run used), Darkwave Disco.
Disclosed differences: individual jobs instead of a Session (pay-as-you-go bills an open Session); backend named instead of least_busy;
the cost Hamiltonian is built directly (0.5 Z_iZ_j per edge, identical to to_ising of the MaxCut QUBO; qiskit-optimization not installed).
FORECAST: Qiskit Aer NoiseModel.from_backend (IBM's record at firing) on the same compiled circuit with the hardware's optimised angles,
200,000 shots (no twirling or decoupling: the standard designer's simulation).
KNOWN-GOOD: the identical code on IBM's local testing mode (AerSimulator.from_backend(ibm_kingston)) reproduces the tutorial: most
likely bitstring cuts 5 (the optimum), 30 COBYLA evaluations (the tutorial reports nfev 30), P(optimal cut) 0.319 against forecast 0.315.
READINGS per chip:
  R1 the tutorial's answer: does the hardware's most likely bitstring give the optimal cut (5)? (YES / NO)
  R2 P(optimal cut), hardware against forecast: 95 % binomial interval of the hardware value against the forecast value -> ABOVE /
     BELOW / CONSISTENT (forecast shot noise at 200,000 shots is below 0.002 and is added to the width).
  R3 mean cut, hardware against forecast, bootstrap over the 10,000 shots, four outcomes.
  R4 described: evaluations, final optimiser cost on hardware, IBM seconds charged.
Use: QIA demonstration (indahl.ai) with every outcome printed; a deviation found here is reported against IBM's own tutorial procedure.
