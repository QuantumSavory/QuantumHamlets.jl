#!/usr/bin/env python3
"""Disposable matching worker with external stage deadlines; writes no audit logs.

Run from any directory: python3 test/run_matching_stress.py [--project=...]
The default project is the package root. The existing experiment environment
can be selected with --project=experiments/bury_performance.
"""
import argparse
import os
from pathlib import Path
import queue
import signal
import subprocess
import sys
import threading
import time

ROOT = Path(__file__).resolve().parent.parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", default=".")
    args = parser.parse_args()
    process = subprocess.Popen(
        ["julia", "--startup-file=no", "--threads=1", f"--project={args.project}",
         str(ROOT / "test" / "matching_stress.jl")],
        cwd=ROOT, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT, text=True, bufsize=1, start_new_session=True,
    )
    messages = queue.Queue()

    def read_output():
        for line in process.stdout:
            messages.put(line.rstrip("\n"))
        messages.put(None)

    threading.Thread(target=read_output, daemon=True).start()

    def until(prefix, seconds):
        deadline = time.monotonic() + seconds
        while True:
            try:
                line = messages.get(timeout=max(0.001, deadline - time.monotonic()))
            except queue.Empty as error:
                raise TimeoutError(f"{prefix!r} exceeded {seconds}s external limit") from error
            if line is None:
                raise RuntimeError("matching worker exited before completing its checks")
            print(line, flush=True)
            if line.startswith(prefix):
                return line

    def stage(kind, n, timeout):
        process.stdin.write(f"{kind} {n}\n")
        process.stdin.flush()
        result = until(f"RESULT\t{kind}\t{n}\t", timeout)
        return float(result.split("\t")[-1])

    try:
        until("READY", 120)
        for kind in ("scalar", "witness"):
            if stage(kind, 1_000, 5) > 1:
                raise RuntimeError("small matching check exceeded 1s; no larger case attempted")
        # Keep a conservative deadline across machines and dependency versions.
        # Never escalate beyond the original stack-overflow reproducer.
        for kind in ("scalar", "witness"):
            stage(kind, 50_000, 60)
        process.stdin.write("QUIT\n")
        process.stdin.flush()
        if process.wait(timeout=5) != 0:
            raise RuntimeError("matching worker failed")
        print("Matching stress regression passed: 50,000 vertices, cost 25,000, valid witness.")
    finally:
        if process.poll() is None:
            os.killpg(process.pid, signal.SIGKILL)
            process.wait()


if __name__ == "__main__":
    try:
        main()
    except (TimeoutError, RuntimeError, subprocess.TimeoutExpired) as error:
        print(f"Matching stress regression failed: {error}", file=sys.stderr)
        raise SystemExit(1)
