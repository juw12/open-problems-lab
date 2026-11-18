"""Extract text from the two OpenAI preprint PDFs to find the exact claimed exponent."""
import os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
PRE = os.path.join(HERE, "openai_preprints")

def extract_pdf_text(path):
    try:
        from pypdf import PdfReader
    except ImportError:
        try:
            from PyPDF2 import PdfReader
        except ImportError:
            print("NO PDF LIBRARY available. pip install pypdf")
            sys.exit(1)
    reader = PdfReader(path)
    text = "\n".join((page.extract_text() or "") for page in reader.pages)
    return text

for name in ("r5t.pdf", "fixed_s.pdf"):
    p = os.path.join(PRE, name)
    print(f"\n{'='*72}\n{name}: {os.path.getsize(p)} bytes\n{'='*72}")
    text = extract_pdf_text(p)
    print(f"  total chars: {len(text)}, pages: {text.count(chr(10))//40}")  # rough
    # Print first 4000 chars (abstract + intro typically)
    print(text[:4000])
    print(f"\n  ... [truncated, {len(text)} total chars] ...")
    # Also print last 2000 chars (often contains the main theorem statement)
    if len(text) > 6000:
        print(f"\n--- TAIL (last 2000 chars) ---\n{text[-2000:]}")
