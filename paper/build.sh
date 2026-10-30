#!/bin/bash
# Builds main.pdf from main.tex with tectonic. The title block carries a fixed date (the date of this version), and SOURCE_DATE_EPOCH pins the PDF's internal timestamps.
# Needs tectonic (https://tectonic-typesetting.github.io). usage: ./build.sh   then   python3 verify_paper.py
cd "$(dirname "$0")" && SOURCE_DATE_EPOCH=1792929600 tectonic main.tex && echo "built main.pdf"
