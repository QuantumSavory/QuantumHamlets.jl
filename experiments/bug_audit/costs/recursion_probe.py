#!/usr/bin/env python3
"""A five-second external stage watchdog, separate from Julia warmup."""
import os
from pathlib import Path
import queue
import signal
import subprocess
import sys
import threading
import time

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
FULL_TARGET=int(sys.argv[1]) if len(sys.argv)>1 else None
assert FULL_TARGET is None or 1 <= FULL_TARGET <= 100000
p=subprocess.Popen(['julia','--startup-file=no','--threads=1',
    '--project=experiments/bury_performance',str(HERE/'recursion_worker.jl')],
    cwd=ROOT,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
    text=True,bufsize=1,start_new_session=True)
lines=queue.Queue()
def read():
    for line in p.stdout: lines.put(line.rstrip('\n'))
    lines.put(None)
threading.Thread(target=read,daemon=True).start()

def kill():
    if p.poll() is None:
        os.killpg(p.pid,signal.SIGKILL)
        p.wait()

with (HERE/('recursion_confirm.log' if FULL_TARGET else 'recursion.log')).open('w') as log:
    def emit(line):
        log.write(line+'\n');log.flush()
        print(line,flush=True)
    def get(deadline):
        try: line=lines.get(timeout=max(.001,deadline-time.monotonic()))
        except queue.Empty:
            kill();raise TimeoutError('external deadline')
        if line is None: raise RuntimeError('worker exited')
        emit(line)
        return line
    def job(kind,n):
        p.stdin.write(f'{kind}\t{n}\n');p.stdin.flush()
        deadline=time.monotonic()+5.0
        while True:
            try: line=get(deadline)
            except TimeoutError:
                emit(f'TIMEOUT\t{kind}\t{n}\tfive_second_limit')
                return None
            if line.startswith('START\t'):
                deadline=time.monotonic()+5.0
            if line.startswith('RESULT\t'): return line.split('\t')
    try:
        deadline=time.monotonic()+120
        while get(deadline)!='READY': pass
        failed=FULL_TARGET
        if FULL_TARGET is None:
            for n in (1000,5000,10000,20000,50000,100000):
                row=job('copied',n)
                if row is None: break
                if row[3]=='StackOverflowError':
                    failed=n;break
                if row[3]!='OK' or float(row[4])>1:
                    emit('STOP\tcopied_probe_slow_or_unexpected_failure');break
        if failed is not None:
            small=job('full',1000)
            if small and small[3]=='OK' and float(small[4])<=1.0:
                job('full',failed)
            else:
                emit('STOP\tno_escalation_after_slow_or_failed_small_full_call')
        else:
            emit('NO_COPIED_STACK_FAILURE_CONFIRMED_WITHIN_LIMIT')
    finally:
        if p.poll() is None:
            p.stdin.write('QUIT\n');p.stdin.flush()
            try:p.wait(timeout=5)
            except subprocess.TimeoutExpired:kill()
