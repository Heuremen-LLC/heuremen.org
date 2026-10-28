"""23 Oct 2026 morning: new record pairs for the restamp check. Saves the records our Hahn jobs saved at firing on 22 Oct (copied from
experiments/2026-10-22-hahn-t2) and IBM's current records for kingston, marrakesh and fez, pulled from the published properties (no job)."""
import json, os, sys, gzip, shutil
HERE = os.path.dirname(os.path.abspath(__file__)); H = os.path.join(HERE, '..', '..', '2026-10-22-hahn-t2')
for src, dst in [('props_hahn_ibm_kingston.json', 'ibm_kingston_2026-10-22_1845.json.gz'), ('props_hahn_ibm_marrakesh.json', 'ibm_marrakesh_2026-10-22_1821.json.gz'),
                 ('props_hahn_ibm_fez_run1.json', 'ibm_fez_2026-10-22_1817.json.gz'), ('props_hahn_ibm_fez_clock1.json', 'ibm_fez_2026-10-22_1955.json.gz'),
                 ('props_hahn_ibm_fez_clock2.json', 'ibm_fez_2026-10-22_2149.json.gz')]:
    json.dump(json.load(open(os.path.join(H, src))), gzip.open(os.path.join(HERE, dst), 'wt'))
sys.path.insert(0, os.path.expanduser('~/bones-forge/arxiv_one_readout/hahn_t2')); import hahn as R
s = R.service()
for chip in ('ibm_kingston', 'ibm_marrakesh', 'ibm_fez'):
    p = s.backend(chip).properties(); d = p.last_update_date.strftime('%Y-%m-%d_%H%M')
    json.dump(p.to_dict(), gzip.open(os.path.join(HERE, f'{chip}_{d}.json.gz'), 'wt'), default=str); print('saved', chip, p.last_update_date)
