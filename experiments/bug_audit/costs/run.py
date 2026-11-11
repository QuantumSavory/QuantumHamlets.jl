#!/usr/bin/env python3
"""Bound this small-graph audit externally, including native matching calls."""
import os
from pathlib import Path
import signal
import subprocess

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
with (HERE/'oracles.log').open('w') as log:
    p = subprocess.Popen(['julia','--startup-file=no','--threads=1',
        '--project=experiments/bury_performance',str(HERE/'oracles.jl')],
        cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
    try:
        code=p.wait(timeout=120)
    except subprocess.TimeoutExpired:
        os.killpg(p.pid,signal.SIGKILL)
        p.wait()
        print('Cost oracle audit exceeded 120s and was stopped.')
        raise SystemExit(124)
print('Cost oracle audit exit:',code)
raise SystemExit(code)
