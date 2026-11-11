#!/usr/bin/env python3
"""Serial Julia benchmarks with external per-stage and per-family guards."""
import csv
import json
import os
from pathlib import Path
import queue
import signal
import subprocess
import threading
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
STAGE_LIMIT = 8.0
STARTUP_LIMIT = 120.0
TOTAL_LIMIT = 600.0


class Worker:
    def __init__(self, log):
        self.log = log
        self.lines = queue.Queue()
        self.process = subprocess.Popen(
            ["julia", "--startup-file=no", "--threads=1",
             "--project=experiments/bury_performance", str(HERE / "worker.jl")],
            cwd=ROOT, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, text=True, bufsize=1, start_new_session=True,
            env=dict(os.environ, OPENBLAS_NUM_THREADS="1"))
        def read():
            for line in self.process.stdout:
                self.lines.put(line.rstrip("\n"))
            self.lines.put(None)
        threading.Thread(target=read, daemon=True).start()
        deadline = time.monotonic() + STARTUP_LIMIT
        while True:
            line = self.get(deadline)
            if line == "READY":
                break

    def get(self, deadline):
        try:
            line = self.lines.get(timeout=max(0.001, deadline-time.monotonic()))
        except queue.Empty:
            self.close(kill=True)
            raise TimeoutError("worker stage exceeded deadline")
        if line is None:
            raise RuntimeError("Julia worker exited unexpectedly")
        self.log.write(line + "\n")
        self.log.flush()
        return line

    def run(self, job):
        self.process.stdin.write("\t".join(str(job[k]) for k in ("job", "topology", "kind", "target", "seed")) + "\n")
        self.process.stdin.flush()
        deadline = time.monotonic() + STAGE_LIMIT
        phase = "dispatch"
        while True:
            try:
                line = self.get(deadline)
            except TimeoutError:
                return "timeout", phase
            fields = line.split("\t")
            if fields[0] == "START":
                phase = "/".join(fields[2:])
                deadline = time.monotonic() + STAGE_LIMIT
            elif fields[0] in ("DONE", "SLOW", "ERROR"):
                return fields[0].lower(), "/".join(fields[2:])

    def close(self, kill=False):
        if self.process.poll() is None:
            if kill:
                os.killpg(self.process.pid, signal.SIGKILL)
            else:
                self.process.stdin.write("QUIT\n")
                self.process.stdin.flush()
            try:
                self.process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                os.killpg(self.process.pid, signal.SIGKILL)
                self.process.wait()


def main():
    with (HERE / "inventory/inventory.tsv").open() as f:
        graphs = list(csv.DictReader(f, delimiter="\t"))
    # Parameter order avoids lexicographic 10-before-3 escalation. Each circuit's
    # initial and reduced forms are visited before the next circuit size.
    graphs.sort(key=lambda r: (r["family"], int(r["parameter"] or 0), r["directory"]))
    representative = {"barenco_tof_3", "barenco_tof_10", "gf2^6_mult", "qcla_com_7", "k1", "k3", "k5"}
    jobs = []
    for r in graphs:
        n = int(r["n"])
        modes = [("ordinary", n//2), ("seeded", n//2), ("automatic", n//2)]
        if r["graph"] in representative:
            modes += [(kind, target) for target in (n//4, 3*n//4) for kind in ("ordinary", "seeded")]
        for kind, target in modes:
            jobs.append(dict(job=f'{r["directory"]}__{r["graph"]}__{kind}__{target}',
                             graph=r["graph"], directory=r["directory"], family=r["family"],
                             parameter=int(r["parameter"] or 0), n=n,
                             topology=str(ROOT/r["topology_path"]), kind=kind, target=target, seed=1))
    (HERE / "plan.json").write_text(json.dumps(jobs, indent=2)+"\n")
    blocked = {}
    failed_graphs = set()
    worker = None
    completed = set()
    start = time.monotonic()
    with (HERE / "worker.log").open("w") as log, (HERE / "status.tsv").open("w") as status:
        writer = csv.DictWriter(status, delimiter="\t", fieldnames=["job","status","detail","wall_seconds"])
        writer.writeheader()
        try:
            worker = Worker(log)
            start = time.monotonic()
            for job in jobs:
                graphkey = (job["directory"], job["graph"])
                reason = None
                if graphkey in failed_graphs:
                    reason = "earlier_workload_on_this_graph_was_slow_or_failed"
                if job["family"] in blocked and job["parameter"] > blocked[job["family"]]:
                    reason = "larger_family_member_after_slow_or_failed_call"
                if time.monotonic()-start > TOTAL_LIMIT:
                    reason = "total_experiment_budget"
                tick = time.monotonic()
                if reason:
                    result, detail = "skipped", reason
                else:
                    for name in ("results.tsv", "samples.tsv", "validation.tsv"):
                        (HERE/"jobs"/job["job"]/name).unlink(missing_ok=True)
                    if worker is None:
                        worker = Worker(log)
                    result, detail = worker.run(job)
                    if result != "done":
                        failed_graphs.add(graphkey)
                        blocked[job["family"]] = min(blocked.get(job["family"], job["parameter"]), job["parameter"])
                    if result == "timeout":
                        worker.close(kill=True)
                        worker = None
                writer.writerow(dict(job=job["job"], status=result, detail=detail, wall_seconds=time.monotonic()-tick))
                status.flush()
                print(job["job"], result, detail, flush=True)
                if result == "done":
                    completed.add(job["job"])
                if result == "error":
                    raise RuntimeError(f'Correctness/worker error for {job["job"]}: {detail}')
        finally:
            if worker is not None:
                worker.close()
    for name in ("results.tsv", "samples.tsv"):
        paths = [HERE/"jobs"/j["job"]/name for j in jobs if j["job"] in completed]
        with (HERE/name).open("w") as out:
            first = True
            for path in paths:
                if not path.exists():
                    continue
                lines = path.read_text().splitlines()
                out.write("\n".join(lines if first else lines[1:])+"\n")
                first = False
    print("GRAPHML_BENCHMARKS_COMPLETE", flush=True)


if __name__ == "__main__":
    main()
