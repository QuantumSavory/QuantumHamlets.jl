#!/usr/bin/env python3
"""Generate GraphML METIS/BURY plots with external time limits and fresh output.

Run from the package root:
    python3 resources/run_graphml_plots.py

Uses the installed parent JuliaProjects environment without resolving packages.
Partial plots, tables, and output.log are retained if a stage fails or times out.
"""
import argparse
from datetime import datetime, timezone
import os
from pathlib import Path
import queue
import shlex
import signal
import subprocess
import sys
import threading
import time


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
STARTUP_SECONDS = 180
STAGE_SECONDS = 90
TOTAL_SECONDS = 900


def objectives(value):
    names = [item.strip() for item in value.split(',')]
    if not names or any(name not in ('matching', 'rank') for name in names):
        raise argparse.ArgumentTypeError('objectives must be matching, rank, or matching,rank')
    if len(set(names)) != len(names):
        raise argparse.ArgumentTypeError('objectives must not contain duplicates')
    return ','.join(names)


def parts(value):
    try:
        counts = [int(item.strip()) for item in value.split(',')]
    except ValueError as error:
        raise argparse.ArgumentTypeError('parts must be comma-separated integers') from error
    if not counts or any(count < 2 for count in counts):
        raise argparse.ArgumentTypeError('each part count must be at least two')
    if len(set(counts)) != len(counts):
        raise argparse.ArgumentTypeError('part counts must not contain duplicates')
    return ','.join(str(count) for count in counts)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--output-dir', type=Path, help='fresh directory; existing paths are rejected')
    parser.add_argument('--objectives', type=objectives, default='matching,rank',
                        help='comma-separated matching and/or rank (default: matching,rank)')
    parser.add_argument('--backend', choices=('heap', 'scan', 'legacy'), default='heap')
    parser.add_argument('--parts', type=parts, default='2,3,4,5,6,7,8',
                        help='comma-separated part counts (default: 2,3,4,5,6,7,8)')
    args = parser.parse_args()
    output = (args.output_dir or HERE / 'figs' / 'graphml_comparison' /
              ('run_' + datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ'))).resolve()
    output.mkdir(parents=True, exist_ok=False)
    command = ['julia', '--startup-file=no', '--threads=1', '--project=..',
               str(HERE / 'graphml_plots.jl'), '--output-dir', str(output),
               '--objectives', args.objectives, '--backend', args.backend, '--parts', args.parts]
    print(f'Output: {output}', flush=True)
    with (output / 'output.log').open('x', encoding='utf-8') as log:
        print('COMMAND\t' + shlex.join(command), file=log, flush=True)
        print(f'LIMITS\tstartup={STARTUP_SECONDS}s\tstage={STAGE_SECONDS}s\ttotal={TOTAL_SECONDS}s',
              file=log, flush=True)
        process = None
        try:
            started = time.monotonic()
            process = subprocess.Popen(command, cwd=ROOT, stdout=subprocess.PIPE,
                                       stderr=subprocess.STDOUT, text=True, encoding='utf-8',
                                       errors='replace', bufsize=1, start_new_session=True)
            lines = queue.Queue()

            def reader():
                for line in process.stdout:
                    lines.put(line.rstrip('\n'))
                lines.put(None)

            threading.Thread(target=reader, daemon=True).start()
            deadline = started + STARTUP_SECONDS
            total_deadline = started + TOTAL_SECONDS
            stage = 'package load'
            complete = False
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
                print(line, flush=True)
                if line == 'READY' or line.startswith('STAGE\t'):
                    stage = line
                    deadline = time.monotonic() + STAGE_SECONDS
                elif line == 'COMPLETE':
                    complete = True
            code = process.wait(timeout=5)
            if code != 0 or not complete:
                raise RuntimeError(f'Julia worker failed (exit {code}, COMPLETE={complete})')
            print(f'ELAPSED_SECONDS\t{time.monotonic() - started:.3f}', file=log, flush=True)
        except (OSError, TimeoutError, RuntimeError, subprocess.TimeoutExpired, KeyboardInterrupt) as error:
            print(f'RUNNER_FAILURE\t{type(error).__name__}: {error}', file=log, flush=True)
            print(f'Partial output retained: {output}', file=sys.stderr, flush=True)
            raise
        finally:
            if process is not None and process.poll() is None:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
    print(f'GraphML plot generation completed: {output}', flush=True)


if __name__ == '__main__':
    try:
        main()
    except (OSError, TimeoutError, RuntimeError, subprocess.TimeoutExpired) as error:
        print(error, file=sys.stderr)
        raise SystemExit(1)
    except KeyboardInterrupt:
        print('Interrupted; partial output was preserved.', file=sys.stderr)
        raise SystemExit(130)
