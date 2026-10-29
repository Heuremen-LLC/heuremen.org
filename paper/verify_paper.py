#!/usr/bin/env python3
"""verify_paper.py: the rebuilt main.pdf must carry the same text as the published paper (https://heuremen.org/a-readout-is-not-an-idle.pdf), ignoring spacing and ligature glyphs, and the same page count.
usage: python3 verify_paper.py [main.pdf] [published.pdf or URL]. Needs pypdf."""
import os, re, sys, tempfile, unicodedata, urllib.request
import pypdf
def txt(p): r = pypdf.PdfReader(p); return len(r.pages), re.sub(r'\s+', '', unicodedata.normalize('NFKC', ''.join((pg.extract_text() or '') for pg in r.pages)))
a = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), 'main.pdf'); pub = sys.argv[2] if len(sys.argv) > 2 else 'https://heuremen.org/a-readout-is-not-an-idle.pdf'
if pub.startswith('http'): p = os.path.join(tempfile.mkdtemp(), 'pub.pdf'); urllib.request.urlretrieve(pub, p); pub = p
n1, t1 = txt(a); n2, t2 = txt(pub); same = n1 == n2 and t1 == t2; print('pages', n1, n2, 'text', len(t1), len(t2), 'IDENTICAL' if same else 'DIFFERENT')
if not same:
    for i, (x, y) in enumerate(zip(t1, t2)):
        if x != y: print('first difference at', i, repr(t1[max(0, i - 40):i + 40]), '|', repr(t2[max(0, i - 40):i + 40])); break
sys.exit(0 if same else 1)
