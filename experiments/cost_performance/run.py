#!/usr/bin/env python3
"""Bound the cost comparison externally and preserve every completed run.

Usage: python3 experiments/cost_performance/run.py [--output-dir PATH]
"""
import argparse
from datetime import datetime, timezone
import os
from pathlib import Path
import queue
import signal
import subprocess
import sys
import threading
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path)
    parser.add_argument('--operation', action='append', choices=(
        'edge_cut_legacy', 'cover_legacy', 'rank_legacy', 'matching_legacy',
        'edge_cut_scalar', 'cover_scalar', 'rank_scalar', 'matching_scalar'),
        help='measure only this operation; repeat to select several (default: all)')
    args = parser.parse_args()
    output = args.output_dir or HERE / ('run_' + datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ'))
    output.mkdir(parents=True, exist_ok=False)
    command = ['julia', '--startup-file=no', '--threads=1',
               '--project=experiments/bury_performance', str(HERE / 'benchmark.jl'), *(args.operation or [])]
    process = subprocess.Popen(command, cwd=ROOT, stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT, text=True, bufsize=1,
                               start_new_session=True)
    lines = queue.Queue()

    def reader():
        for line in process.stdout:
            lines.put(line.rstrip('\n'))
        lines.put(None)

    threading.Thread(target=reader, daemon=True).start()
    deadline = time.monotonic() + 180
    total_deadline = time.monotonic() + 480
    stage = 'package load'
    complete = False
    try:
        with (output / 'output.log').open('x') as log, (output / 'results.tsv').open('x') as results:
            while True:
                remaining = min(deadline, total_deadline) - time.monotonic()
                if remaining <= 0:
                    raise TimeoutError(f'{stage} exceeded its external deadline')
                try:
                    line = lines.get(timeout=remaining)
                except queue.Empty as error:
                    raise TimeoutError(f'{stage} exceeded its external deadline') from error
                if line is None:
                    break
                print(line, file=log, flush=True)
                if line.startswith('HEADER\t') or line.startswith('RESULT\t'):
                    print(line.split('\t', 1)[1], file=results, flush=True)
                if line == 'READY' or line.startswith('STAGE\t'):
                    stage = line
                    deadline = time.monotonic() + 60
                    print(line, flush=True)
                elif line.startswith('RESULT\t'):
                    print(line, flush=True)
                elif line == 'COMPLETE':
                    complete = True
            code = process.wait(timeout=5)
            if code != 0 or not complete:
                raise RuntimeError(f'worker failed (exit {code}); see {output / "output.log"}')
        print(f'Cost comparison passed; output: {output}')
    finally:
        if process.poll() is None:
            os.killpg(process.pid, signal.SIGKILL)
            process.wait()


if __name__ == '__main__':
    try:
        main()
    except (TimeoutError, RuntimeError, subprocess.TimeoutExpired) as error:
        print(error, file=sys.stderr)
        raise SystemExit(1)
