"""
Capture and print computational environment information.

Usage:
    PYTHONPATH=src python scripts/capture_session_info.py
    PYTHONPATH=src python scripts/capture_session_info.py --json logs/session_info.json

Writes a machine-readable JSON manifest and a human-readable text report
to stdout / the specified paths.  Used for the reproducibility appendix and
to verify environments before / after a pipeline run.
"""

from __future__ import annotations

import argparse
import importlib.metadata as importlib_metadata
import json
import platform
import sys
import time
from pathlib import Path

KEY_PACKAGES = [
    "meningeal-extension",
    "numpy",
    "pandas",
    "scipy",
    "statsmodels",
    "scikit-learn",
    "matplotlib",
    "seaborn",
    "plotly",
    "GEOparse",
    "mygene",
    "biopython",
    "lifelines",
    "tqdm",
    "pytest",
]


def get_versions(packages: list[str]) -> dict[str, str | None]:
    out: dict[str, str | None] = {}
    for pkg in packages:
        try:
            out[pkg] = importlib_metadata.version(pkg)
        except importlib_metadata.PackageNotFoundError:
            try:
                mod = importlib_metadata.metadata(pkg.replace("-", "_"))
                out[pkg] = mod.get("Version", None)
            except Exception:
                out[pkg] = None
    return out


def collect_info() -> dict:
    info = {
        "timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "python": {
            "version": sys.version,
            "executable": sys.executable,
            "implementation": platform.python_implementation(),
        },
        "platform": {
            "system": platform.system(),
            "release": platform.release(),
            "version": platform.version(),
            "machine": platform.machine(),
            "processor": platform.processor(),
        },
        "packages": get_versions(KEY_PACKAGES),
    }
    return info


def format_text(info: dict) -> str:
    lines = []
    lines.append("=" * 64)
    lines.append("COMPUTATIONAL ENVIRONMENT — SESSION INFO")
    lines.append("=" * 64)
    lines.append(f"Timestamp (UTC) : {info['timestamp_utc']}")
    lines.append(f"Python          : {info['python']['version'].splitlines()[0]}")
    lines.append(f"  executable    : {info['python']['executable']}")
    lines.append(f"  implementation: {info['python']['implementation']}")
    lines.append(
        f"Platform        : {info['platform']['system']} "
        f"{info['platform']['release']} ({info['platform']['machine']})"
    )
    lines.append("")
    lines.append("Key package versions:")
    lines.append("-" * 40)
    for pkg, ver in sorted(info["packages"].items()):
        marker = "  [MISSING]" if ver is None else ""
        lines.append(f"  {pkg:<24s} : {ver or 'NOT INSTALLED'}{marker}")
    lines.append("")
    lines.append("Full pip freeze: see requirements-freeze.txt at repo root.")
    lines.append("=" * 64)
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--json",
        type=Path,
        default=None,
        help="Write machine-readable JSON to this path.",
    )
    parser.add_argument(
        "--text",
        type=Path,
        default=None,
        help="Write human-readable text report to this path.",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Do not print the text report to stdout.",
    )
    args = parser.parse_args()

    info = collect_info()

    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(info, indent=2, sort_keys=True), encoding="utf-8")
        if not args.quiet:
            print(f"[session_info] JSON manifest written to: {args.json}")

    report = format_text(info)
    if args.text:
        args.text.parent.mkdir(parents=True, exist_ok=True)
        args.text.write_text(report + "\n", encoding="utf-8")
        if not args.quiet:
            print(f"[session_info] Text report written to: {args.text}")

    if not args.quiet:
        print()
        print(report)


if __name__ == "__main__":
    main()
