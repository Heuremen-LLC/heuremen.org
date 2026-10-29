#!/bin/bash
# Builds main.pdf from main.tex with tectonic. The title block uses \today, so the build date is pinned to the date the paper was published (25 Oct 2026, 12:00 UTC) with SOURCE_DATE_EPOCH.
# Needs tectonic (https://tectonic-typesetting.github.io). usage: ./build.sh   then   python3 verify_paper.py
cd "$(dirname "$0")" && SOURCE_DATE_EPOCH=1792929600 tectonic main.tex && echo "built main.pdf"
