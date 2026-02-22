#!/usr/bin/env python3
import fitz
import sys

if len(sys.argv) < 2:
    print("Usage: python3 read_pdf2.py <pdf_file> [start_page] [end_page]")
    sys.exit(1)

pdf_file = sys.argv[1]
start_page = int(sys.argv[2]) if len(sys.argv) > 2 else 0
end_page = int(sys.argv[3]) if len(sys.argv) > 3 else start_page + 10

doc = fitz.open(pdf_file)
print(f"=== {pdf_file} ===")
print(f"Total pages: {len(doc)}")
print(f"Reading pages {start_page+1} to {min(end_page, len(doc))}")
print()

for i in range(start_page, min(end_page, len(doc))):
    page = doc[i]
    print(f"--- Page {i+1} ---")
    print(page.get_text()[:8000])
    print()
