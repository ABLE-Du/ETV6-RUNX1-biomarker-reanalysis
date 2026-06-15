from __future__ import annotations

import csv
import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT / "RELEASE_SHA256_MANIFEST.csv"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    files = sorted(
        path for path in ROOT.rglob("*")
        if path.is_file() and path != OUTPUT and "__pycache__" not in path.parts
    )
    with OUTPUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["relative_path", "size_bytes", "sha256"])
        for path in files:
            writer.writerow([path.relative_to(ROOT).as_posix(), path.stat().st_size, sha256(path)])


if __name__ == "__main__":
    main()
