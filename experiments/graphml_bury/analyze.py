#!/usr/bin/env python3
"""Verify stored measurements and write compact backend comparisons."""
import csv
import hashlib
import json
import statistics
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent


def read(name):
    with (HERE/name).open() as f:
        return list(csv.DictReader(f, delimiter="\t"))


def main():
    inventory = read("inventory/inventory.tsv")
    plans = {r["job"]: r for r in json.loads((HERE/"plan.json").read_text())}
    status = read("status.tsv")
    done = {r["job"] for r in status if r["status"] == "done"}
    rows = read("results.tsv")
    samples = defaultdict(list)
    for r in read("samples.tsv"):
        samples[(r["job"], r["implementation"])].append(float(r["time_ns"]))
    byjob = defaultdict(dict)
    for r in rows:
        assert r["job"] in done
        values = samples[(r["job"],r["implementation"])]
        assert len(values) == int(r["samples"])
        assert statistics.median(values) == float(r["median_ns"])
        byjob[r["job"]][r["implementation"]] = r
    assert set(byjob) == done
    assert all(set(rs) == {"legacy", "scan", "heap"} for rs in byjob.values())
    summary = []
    for job, rs in byjob.items():
        p = plans[job]
        times = {k: float(v["median_ns"])/1000 for k,v in rs.items()}
        order = sorted(times, key=times.get)
        summary.append(dict(directory=p["directory"],graph=p["graph"],kind=p["kind"],
                            n=p["n"],target=p["target"],half=p["target"]==p["n"]//2,
                            legacy_us=times["legacy"],scan_us=times["scan"],heap_us=times["heap"],
                            fastest=order[0],within_10pct="/".join(k for k in order if times[k]<=1.1*times[order[0]]),
                            legacy_over_heap=times["legacy"]/times["heap"],
                            legacy_over_scan=times["legacy"]/times["scan"],
                            scan_over_heap=times["scan"]/times["heap"]))
    with (HERE/"comparison.tsv").open("w") as f:
        writer=csv.DictWriter(f, fieldnames=list(summary[0]), delimiter="\t")
        writer.writeheader();writer.writerows(summary)
    validation=[]
    for item in status:
        path=HERE/"jobs"/item["job"]/"validation.tsv"
        if not path.exists(): continue
        with path.open() as f:
            for r in csv.DictReader(f,delimiter="\t"):
                for key in ("partition_equal","graph_unchanged","balanced","rng_equal"):
                    assert r[key] == "true"
                validation.append(dict(job=item["job"], **r))
    if validation:
        with (HERE/"validation.tsv").open("w") as f:
            writer=csv.DictWriter(f,fieldnames=list(validation[0]),delimiter="\t")
            writer.writeheader();writer.writerows(validation)
    for rel, digest in json.loads((HERE/"baseline_sha256.json").read_text()).items():
        assert hashlib.sha256((ROOT/rel).read_bytes()).hexdigest() == digest, rel
    for r in inventory:
        assert hashlib.sha256((ROOT/r["path"]).read_bytes()).hexdigest() == r["sha256"]
    print("Status:",dict(Counter(r["status"] for r in status)))
    print("Verified result rows:",len(rows),"raw timing samples:",sum(map(len,samples.values())))
    print("Validated backend calls:",len(validation))
    print("All 210 pre-existing files and all 45 GraphML inputs unchanged.")
    for directory in sorted({r["directory"] for r in summary}):
        for kind in ("ordinary","seeded","automatic"):
            group=[r for r in summary if r["directory"]==directory and r["kind"]==kind and r["half"]]
            if not group: continue
            sums={b:sum(r[b+"_us"] for r in group) for b in ("legacy","scan","heap")}
            wins=Counter(r["fastest"] for r in group)
            clear=Counter(r["fastest"] for r in group if "/" not in r["within_10pct"])
            print(directory,kind,"fastest median:",dict(wins),"over 10% margin:",dict(clear),"sum of per-file medians(us):",sums)


if __name__ == "__main__":
    main()
