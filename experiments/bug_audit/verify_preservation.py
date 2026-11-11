#!/usr/bin/env python3
import hashlib
import json
from pathlib import Path

here = Path(__file__).resolve().parent
root = here.parent.parent
baseline = json.loads((here / "baseline_sha256.json").read_text())
changed = [name for name, digest in baseline.items()
           if not (root / name).is_file()
           or hashlib.sha256((root / name).read_bytes()).hexdigest() != digest]
if changed:
    raise SystemExit("Pre-existing files changed: " + ", ".join(changed))
print(f"All {len(baseline)} pre-existing files remain byte-for-byte unchanged.")
