"""Verify pinned hashes for cross-compartment source artifacts."""
from __future__ import annotations

import csv
import hashlib
from pathlib import Path
from typing import Iterable


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _rows(manifest_path: Path) -> list[dict[str, str]]:
    if not manifest_path.exists():
        raise FileNotFoundError(
            f"source hash manifest absent: {manifest_path}. "
            "Pin the cross-compartment source before reading it."
        )
    with manifest_path.open(newline="") as fh:
        return list(csv.DictReader(fh))


def _manifest_key(path: Path, root: Path) -> str:
    try:
        return str(Path(path).resolve().relative_to(root.resolve()))
    except ValueError:
        return str(Path(path))


def verify_source_hash(
    path: Path,
    source_label: str,
    manifest_path: Path,
    *,
    root: Path,
) -> str:
    """Return the current SHA-256 after asserting it matches the pinned manifest row."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"cross-compartment source absent: {path}")
    source_path = _manifest_key(path, root)
    hits = [
        row for row in _rows(manifest_path)
        if row.get("source_label") == source_label and row.get("source_path") == source_path
    ]
    if len(hits) != 1:
        raise RuntimeError(
            f"source hash pin missing or duplicated for {source_label} ({source_path}) "
            f"in {manifest_path}"
        )
    current = sha256_file(path)
    expected = hits[0].get("sha256", "")
    if current != expected:
        raise RuntimeError(
            f"source hash mismatch for {source_label}: {source_path}\n"
            f"expected {expected}\nobserved {current}\n"
            "Regenerate or update the consuming stage only after reviewing the source change."
        )
    return current


def verify_source_hashes(
    manifest_path: Path,
    sources: Iterable[tuple[str, Path]],
    *,
    root: Path,
) -> dict[str, str]:
    return {
        label: verify_source_hash(path, label, manifest_path, root=root)
        for label, path in sources
    }
