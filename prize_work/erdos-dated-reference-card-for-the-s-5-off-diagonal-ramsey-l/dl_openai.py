"""Download the two OpenAI 'family 170' preprints (Sep 24 2026) via proxy.
These are the post-cutoff manuscripts that the latent-seminar page says
'determine the logarithmic exponent' for the s>=5 off-diagonal Ramsey line.
Writes to openai_preprints/ and reports sizes.
"""
import urllib.request, urllib.parse
from pathlib import Path

PROXY = "http://localhost:3128"
HERE = Path(__file__).parent
OUT = HERE / "openai_preprints"
OUT.mkdir(exist_ok=True)

COMMIT = "adc7f1241b42e322a6451854ab7e4b4c146bf78a"
FILES = {
    "r5t": f"preprints/The-Sharp-Logarithmic-Exponent-of-r-5-t-September-24-2026/paper.pdf",
    "fixed_s": f"preprints/Sharp-Logarithmic-Exponents-for-Fixed-Off-Diagonal-Ramsey-Numbers-September-24-2026/paper.pdf",
}

def download(key, path):
    url = f"https://raw.githubusercontent.com/openai/math/{COMMIT}/{path}"
    req = urllib.request.Request(url, headers={"User-Agent": "millennium-lab/1.0"})
    opener = urllib.request.build_opener(urllib.request.ProxyHandler({"http": PROXY, "https": PROXY}))
    dest = OUT / f"{key}.pdf"
    with opener.open(req, timeout=120) as r:
        data = r.read()
    dest.write_bytes(data)
    return dest, len(data)

for key, p in FILES.items():
    try:
        dest, n = download(key, p)
        print(f"{key}: {dest.name} bytes={n}  ({n/1024:.1f} KB)")
    except Exception as ex:
        print(f"{key}: FAILED {type(ex).__name__}: {ex}")
print("DONE")
