"""Verify that every file present before this experiment is byte-for-byte unchanged."""
from pathlib import Path
import hashlib
import json

experiment = Path(__file__).resolve().parent
root = experiment.parent.parent
baseline = json.loads((experiment / "baseline_sha256.json").read_text())
changed = [
    name for name, digest in baseline.items()
    if not (root / name).is_file()
    or hashlib.sha256((root / name).read_bytes()).hexdigest() != digest
]
if changed:
    raise SystemExit("Changed or missing original files:\n" + "\n".join(changed))
print(f"All {len(baseline)} pre-existing files are byte-for-byte unchanged.")
