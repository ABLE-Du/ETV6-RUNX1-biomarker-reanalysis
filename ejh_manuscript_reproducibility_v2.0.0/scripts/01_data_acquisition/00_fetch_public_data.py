#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""00_fetch_public_data.py -- Round-1 acquisition of the open datasets.

Downloads (public, no credentials, no login):
  1. Oksa et al. Leukemia 2025 supplementary package (Europe PMC)  -> data/oksa_supp/
  2. GEO GSE227832 raw read counts + the four per-platform series matrices
     (metadata only, a few kB each) and the GSE228632 counts       -> data/geo/

Transport note: the GEO FTP endpoint throttles a single connection to roughly
35 KB/s on this machine, so every file is fetched with N parallel byte-range
requests and the parts are concatenated.  Parts are assembled only after all of
them arrive, and the concatenation is validated against the server's
Content-Length and with a gzip integrity test.  A partial download can never be
mistaken for a complete one.

Every artefact is recorded in data/PROVENANCE.tsv with url, bytes, sha256,
integrity status and retrieval time.  Existing files are re-verified rather than
re-downloaded, so the script is idempotent.

No analysis is performed here -- Round 1 is acquisition + audit only.
"""
from __future__ import annotations

import concurrent.futures as cf
import datetime as dt
import gzip
import hashlib
import json
import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
OKSA = DATA / "oksa_supp"
GEO = DATA / "geo"
UA = {"User-Agent": "Mozilla/5.0 (ETV6-RUNX1 adverse-response project)"}
N_PARTS = 8

GEO_FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series"
SOURCES = [
    (OKSA, "PMC12380598_SupplementaryFiles.zip",
     "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12380598/supplementaryFiles", True, False),
    (GEO, "GSE227832_RNAseq_read_counts.txt.gz",
     f"{GEO_FTP}/GSE227nnn/GSE227832/suppl/GSE227832_RNAseq_read_counts.txt.gz", True, True),
    (GEO, "GSE228632_RNAseq_read_counts.txt.gz",
     f"{GEO_FTP}/GSE228nnn/GSE228632/suppl/GSE228632_RNAseq_read_counts.txt.gz", True, True),
    *[(GEO, f"GSE227832-{g}_series_matrix.txt.gz",
       f"{GEO_FTP}/GSE227nnn/GSE227832/matrix/GSE227832-{g}_series_matrix.txt.gz", True, True)
      for g in ("GPL11154", "GPL15520", "GPL16791", "GPL24676")],
    (GEO, "GSE228632_series_matrix.txt.gz",
     f"{GEO_FTP}/GSE228nnn/GSE228632/matrix/GSE228632_series_matrix.txt.gz", True, True),
]

# supplementary members we actually need (the package also ships a 19 MB PDF)
KEEP_MEMBERS = {
    "41375_2025_2683_MOESM2_ESM.xlsx": "Oksa2025_supplementary_tables.xlsx",
    "41375_2025_2683_MOESM1_ESM.pdf": "Oksa2025_supplementary_information.pdf",
}


def say(*a):
    print(*a, flush=True)


def head_len(url: str, timeout: int = 60) -> int:
    req = urllib.request.Request(url, headers=UA, method="HEAD")
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return int(r.headers.get("Content-Length", 0))


def get_range(url: str, start: int, end: int, timeout: int = 900) -> bytes:
    h = dict(UA)
    h["Range"] = f"bytes={start}-{end}"
    req = urllib.request.Request(url, headers=h)
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()


def fetch_parallel(url: str, dest: Path, n_parts: int = N_PARTS) -> tuple[int, int]:
    """Download url -> dest with parallel range requests.  Returns (got, want)."""
    want = head_len(url)
    if want <= 0:
        raise RuntimeError("server did not report Content-Length")
    chunk = (want + n_parts - 1) // n_parts
    spans = []
    for i in range(n_parts):
        s = i * chunk
        e = min(s + chunk - 1, want - 1)
        if s <= e:
            spans.append((i, s, e))

    parts: dict[int, bytes] = {}
    with cf.ThreadPoolExecutor(max_workers=n_parts) as ex:
        futs = {ex.submit(get_range, url, s, e): i for i, s, e in spans}
        for fut in cf.as_completed(futs):
            parts[futs[fut]] = fut.result()
    blob = b"".join(parts[i] for i, _, _ in spans)
    dest.write_bytes(blob)
    return len(blob), want


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def gzip_ok(p: Path) -> bool:
    if not p.name.endswith(".gz"):
        return True
    try:
        with gzip.open(p, "rb") as fh:
            while fh.read(1 << 20):
                pass
        return True
    except Exception:  # noqa: BLE001
        return False


def main() -> int:
    OKSA.mkdir(parents=True, exist_ok=True)
    GEO.mkdir(parents=True, exist_ok=True)

    prov: list[dict] = []
    for subdir, name, url, required, is_gz in SOURCES:
        dest = subdir / name
        status = "present"
        if not (dest.exists() and dest.stat().st_size > 0):
            try:
                got, want = fetch_parallel(url, dest)
                say(f"[ok]   {name}: {got:,} B (want {want:,})")
                status = "downloaded" if got == want else "incomplete"
                if got != want:
                    if required:
                        say(f"[FAIL] {name}: size mismatch {got} != {want}")
                        return 1
                    dest.unlink(missing_ok=True)
                    status = "error: size mismatch"
            except Exception as e:  # noqa: BLE001
                say(f"[{'FAIL' if required else 'warn'}] {name}: {e}")
                prov.append(dict(file=str(dest.relative_to(ROOT)), url=url, bytes=0,
                                 sha256="", integrity=f"error: {type(e).__name__}",
                                 status="error"))
                if required:
                    return 1
                continue
        else:
            say(f"[skip] {name} already present ({dest.stat().st_size:,} B)")

        integ = "gzip_ok" if gzip_ok(dest) else "GZIP_FAIL"
        prov.append(dict(file=str(dest.relative_to(ROOT)), url=url,
                         bytes=dest.stat().st_size, sha256=sha256_file(dest),
                         integrity=integ, status=status))

    # ---- unpack the Oksa supplementary members we need -----------------------
    zpath = OKSA / "PMC12380598_SupplementaryFiles.zip"
    if zpath.exists():
        import zipfile
        with zipfile.ZipFile(zpath) as z:
            members = z.namelist()
            say(f"[zip] {len(members)} members")
            for m in members:
                if m in KEEP_MEMBERS:
                    out = OKSA / KEEP_MEMBERS[m]
                    out.write_bytes(z.read(m))
                    say(f"       -> {out.name} ({out.stat().st_size:,} B)")
            (OKSA / "MEMBER_LIST.txt").write_text("\n".join(members) + "\n", encoding="utf-8")

    # ---- verification of the two files the audits depend on ------------------
    checks = []
    x = OKSA / "Oksa2025_supplementary_tables.xlsx"
    checks.append(("oksa supplementary tables xlsx present", x.exists() and x.stat().st_size > 1_000_000))
    c = GEO / "GSE227832_RNAseq_read_counts.txt.gz"
    n_samples = None
    if c.exists():
        with gzip.open(c, "rt", encoding="utf-8", errors="replace") as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
        n_samples = len(hdr) - 1
        checks.append(("GSE227832 counts are gene x sample",
                       n_samples > 100 and hdr[0].lower().lstrip("\ufeff").startswith("gene")))
    checks.append(("no gzip integrity failure in provenance",
                   all(p["integrity"] != "GZIP_FAIL" for p in prov)))

    ts = dt.datetime.now().astimezone().isoformat(timespec="seconds")
    with (DATA / "PROVENANCE.tsv").open("w", encoding="utf-8") as fh:
        fh.write("file\turl\tbytes\tsha256\tintegrity\tstatus\tretrieved_at\n")
        for p in prov:
            fh.write(f"{p['file']}\t{p['url']}\t{p['bytes']}\t{p['sha256']}\t"
                     f"{p['integrity']}\t{p['status']}\t{ts}\n")
    say(f"[write] data/PROVENANCE.tsv ({len(prov)} rows)")

    manifest = {"retrieved_at": ts, "n_files": len(prov),
                "GSE227832_n_samples": n_samples,
                "checks": [{"check": k, "pass": bool(v)} for k, v in checks],
                "all_checks_pass": all(v for _, v in checks)}
    (DATA / "_fetch_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    for k, v in checks:
        say(f"  [{'PASS' if v else 'FAIL'}] {k}")
    say(f"[done] all_checks_pass={manifest['all_checks_pass']}  "
        f"GSE227832 samples={n_samples}")
    return 0 if manifest["all_checks_pass"] else 1


if __name__ == "__main__":
    sys.exit(main())
