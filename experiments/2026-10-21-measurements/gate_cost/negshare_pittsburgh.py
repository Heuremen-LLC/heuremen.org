"""Share of individual qubits with a negative per-gate slope on ibm_pittsburgh (run 121, 2 October), with verify_gate_cost.py's own
slope_stats and sign test, on each of its three rulers. Defined after the data (section 4 of verify_gate_cost.py lists the other chips)."""
import json, os
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'verify_gate_cost.py')).read().split('def main():')[0]
exec(src)
D = json.load(open(os.path.join(HERE, 'data', 'gate_cost.json')))['runs']
for r in D:
    if r['chip'] == 'ibm_pittsburgh':
        for ruler in ('own', 'r90', 'r300'):
            pb = slope_stats(r, ruler)[3]; n = len(pb); neg = sum(1 for x in pb if x < 0)
            print(f"ibm_pittsburgh run {r['run']} ruler {ruler:<4} {n} qubits  {neg} negative = {100*neg/n:.1f} %  (sign-test p = {binom_two_sided(neg, n):.3f})")
