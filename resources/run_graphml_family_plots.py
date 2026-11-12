#!/usr/bin/env python3
"""Generate circuit-family comparisons with fixed METIS seeds and time limits.

Run from the package root:
    python3 resources/run_graphml_family_plots.py

Replot a saved analysis into a fresh directory without recomputing partitions:
    python3 resources/run_graphml_family_plots.py --input-dir PATH

The installed parent JuliaProjects environment is used without resolving packages.
Full worker output and any partial results are preserved in the output directory.
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
STAGE_SECONDS = 120
TOTAL_SECONDS = 1200
PROGRESS_SECONDS = 5


def integer_list(value, *, name, minimum, maximum=None):
    try:
        values = [int(item.strip()) for item in value.split(',')]
    except ValueError as error:
        raise argparse.ArgumentTypeError(f'{name} must be comma-separated integers') from error
    if not values or any(item < minimum for item in values):
        raise argparse.ArgumentTypeError(f'each {name[:-1]} must be at least {minimum}')
    if maximum is not None and any(item > maximum for item in values):
        raise argparse.ArgumentTypeError(f'each {name[:-1]} must be at most {maximum}')
    if len(set(values)) != len(values):
        raise argparse.ArgumentTypeError(f'{name} must not contain duplicates')
    return ','.join(str(item) for item in values)


def seeds(value):
    # The installed METIS wrapper stores its native seed in an Int32 option.
    return integer_list(value, name='seeds', minimum=0, maximum=2**31 - 1)


def parts(value):
    return integer_list(value, name='parts', minimum=2)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--output-dir', type=Path,
                        help='fresh directory; existing paths are rejected')
    parser.add_argument('--input-dir', type=Path,
                        help='existing saved-analysis directory to replot without recomputing')
    parser.add_argument('--seeds', type=seeds, default='0,1,2,3,4,5,6,7,8,9',
                        help='fixed comma-separated METIS seeds (default: 0 through 9)')
    parser.add_argument('--parts', type=parts, default='2,3,4,5,6,7,8',
                        help='comma-separated part counts (default: 2 through 8)')
    parser.add_argument('--backend', choices=('heap', 'scan', 'legacy'), default='heap')
    args = parser.parse_args(argv)
    input_dir = args.input_dir.resolve() if args.input_dir is not None else None
    if input_dir is not None and not input_dir.is_dir():
        parser.error('--input-dir must be an existing directory')
    output = (args.output_dir or HERE / 'figs' / 'graphml_families' /
              ('run_' + datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ'))).resolve()
    output.mkdir(parents=True, exist_ok=False)
    command = ['julia', '--startup-file=no', '--threads=1', '--project=..',
               str(HERE / 'graphml_family_plots.jl'), '--output-dir', str(output),
               '--seeds', args.seeds, '--parts', args.parts, '--backend', args.backend]
    if input_dir is not None:
        command.extend(('--input-dir', str(input_dir)))
    print(f'Output: {output}', flush=True)
    if input_dir is not None:
        print(f'Replot saved analysis: {input_dir}', flush=True)
    with (output / 'output.log').open('x', encoding='utf-8') as log:
        print('COMMAND\t' + shlex.join(command), file=log, flush=True)
        print(f'LIMITS\tstartup={STARTUP_SECONDS}s\tstage={STAGE_SECONDS}s\ttotal={TOTAL_SECONDS}s',
              file=log, flush=True)
        process = None
        try:
            started = time.monotonic()
            process = subprocess.Popen(command, cwd=ROOT, stdin=subprocess.DEVNULL,
                                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                       text=True, encoding='utf-8', errors='replace', bufsize=1,
                                       start_new_session=True)
            lines = queue.Queue()

            def reader():
                try:
                    for line in process.stdout:
                        lines.put(line.rstrip('\n'))
                finally:
                    lines.put(None)

            threading.Thread(target=reader, daemon=True).start()
            deadline = started + STARTUP_SECONDS
            total_deadline = started + TOTAL_SECONDS
            last_progress = started - PROGRESS_SECONDS
            stage = 'package load'
            complete = False
            while True:
                remaining = min(deadline, total_deadline) - time.monotonic()
                if remaining <= 0:
                    cause = 'total run' if time.monotonic() >= total_deadline else stage
                    raise TimeoutError(f'{cause} exceeded its external deadline')
                try:
                    line = lines.get(timeout=remaining)
                except queue.Empty as error:
                    cause = 'total run' if time.monotonic() >= total_deadline else stage
                    raise TimeoutError(f'{cause} exceeded its external deadline') from error
                if line is None:
                    break
                print(line, file=log, flush=True)
                if line == 'READY' or line.startswith('STAGE\t'):
                    stage = line
                    now = time.monotonic()
                    deadline = now + STAGE_SECONDS
                    if line == 'READY' or now - last_progress >= PROGRESS_SECONDS:
                        print(line.replace('\t', '  '), flush=True)
                        last_progress = now
                elif line == 'COMPLETE':
                    complete = True
                    print(line, flush=True)
            code = process.wait(timeout=max(0.01, min(5, total_deadline - time.monotonic())))
            if code != 0 or not complete:
                raise RuntimeError(f'Julia worker failed (exit {code}, COMPLETE={complete}); see output.log')
            print(f'ELAPSED_SECONDS\t{time.monotonic() - started:.3f}', file=log, flush=True)
        except (OSError, TimeoutError, RuntimeError, subprocess.TimeoutExpired, KeyboardInterrupt) as error:
            print(f'RUNNER_FAILURE\t{type(error).__name__}: {error}', file=log, flush=True)
            print(f'Partial output retained: {output}', file=sys.stderr, flush=True)
            raise
        finally:
            if process is not None and process.poll() is None:
                try:
                    os.killpg(process.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                process.wait()
    print(f'GraphML family plots completed: {output}', flush=True)


if __name__ == '__main__':
    try:
        main()
    except (OSError, TimeoutError, RuntimeError, subprocess.TimeoutExpired) as error:
        print(error, file=sys.stderr)
        raise SystemExit(1)
    except KeyboardInterrupt:
        print('Interrupted; partial output was preserved.', file=sys.stderr)
        raise SystemExit(130)
