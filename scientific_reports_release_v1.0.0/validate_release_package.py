from __future__ import annotations

import csv
import hashlib
import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parent
MANIFEST = ROOT / "RELEASE_SHA256_MANIFEST.csv"
EXCLUDED = {MANIFEST}
RISKY_SUFFIXES = {".xlsx", ".xls", ".sav", ".dta", ".db", ".sqlite"}
RISKY_FILENAME_PARTS = {"patient_level", "single_center"}
TEXT_SUFFIXES = {".cff", ".csv", ".json", ".md", ".py", ".txt"}
SENSITIVE_TEXT = re.compile(
    r"TARGET-\d{2}-[A-Z]{6}|C:\\Users\\|E:\\|password\s*=|secret\s*=|api[_-]?key\s*=",
    re.IGNORECASE,
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    files = sorted(
        path for path in ROOT.rglob("*")
        if path.is_file() and path not in EXCLUDED and "__pycache__" not in path.parts
    )
    problems: list[str] = []

    for path in files:
        relative = path.relative_to(ROOT).as_posix()
        lowered = relative.lower()
        if path.suffix.lower() in RISKY_SUFFIXES:
            problems.append(f"disallowed file type: {relative}")
        if any(part in lowered for part in RISKY_FILENAME_PARTS):
            problems.append(f"patient-level filename pattern: {relative}")
        if path.suffix.lower() in TEXT_SUFFIXES and path.name != Path(__file__).name:
            text = path.read_text(encoding="utf-8", errors="replace")
            if SENSITIVE_TEXT.search(text):
                problems.append(f"sensitive text pattern: {relative}")

    json.loads((ROOT / ".zenodo.json").read_text(encoding="utf-8"))

    with MANIFEST.open(encoding="utf-8", newline="") as handle:
        recorded = {row["relative_path"]: row["sha256"] for row in csv.DictReader(handle)}
    current = {path.relative_to(ROOT).as_posix(): sha256(path) for path in files}
    if recorded != current:
        problems.append("SHA-256 manifest does not match current release files")

    if problems:
        raise SystemExit("RELEASE VALIDATION FAILED\n- " + "\n- ".join(problems))
    print(f"RELEASE VALIDATION PASSED: {len(files)} files; no patient-level or sensitive files detected.")


if __name__ == "__main__":
    main()
