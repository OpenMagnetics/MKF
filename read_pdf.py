#!/usr/bin/env python3
import fitz
import sys

if len(sys.argv) < 2:
    print("Usage: python3 read_pdf.py <pdf_file>")
    sys.exit(1)

pdf_file = sys.argv[1]
doc = fitz.open(pdf_file)
print(f"=== {pdf_file} ===")
print(f"Total pages: {len(doc)}")
print()

for i, page in enumerate(doc[:5]):  # First 5 pages
    print(f"--- Page {i+1} ---")
    print(page.get_text()[:5000])
    print()
