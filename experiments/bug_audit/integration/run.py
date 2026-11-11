#!/usr/bin/env python3
"""Bounded probes: no resource sweep is imported or executed."""
import os,signal,subprocess,threading,queue,time
from pathlib import Path
p=Path(__file__).resolve().parent
root=p.parents[2]
q=queue.Queue()
proc=subprocess.Popen(['julia','--startup-file=no','--threads=1','--project=experiments/bury_performance',str(p/'probes.jl')],cwd=root,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,start_new_session=True)
def read():
 for line in proc.stdout:q.put(line)
 q.put(None)
threading.Thread(target=read,daemon=True).start()
deadline=time.monotonic()+120
with (p/'results.log').open('w') as log:
 while True:
  try:line=q.get(timeout=max(.001,deadline-time.monotonic()))
  except queue.Empty:
   os.killpg(proc.pid,signal.SIGKILL);proc.wait()
   log.write('TIMEOUT: active stage stopped\n');print('TIMEOUT');break
  if line is None:
   proc.wait();print('exit',proc.returncode);break
  log.write(line);log.flush()
  if line.startswith('START\t'):
   # First calls may compile solver/plot code; no stage may run indefinitely.
   deadline=time.monotonic()+45
  elif line.startswith(('END\t','READY')):deadline=time.monotonic()+45
