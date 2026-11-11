"""Run the four bounded approximation checks without altering saved evidence."""
from pathlib import Path
import os
import signal
import subprocess
import sys


def main():
    here = Path(__file__).resolve().parent
    root = here.parent.parent
    destination = here / "rerun_logs"
    destination.mkdir(exist_ok=True)
    scripts = (
        "counterexamples/path_hubs.jl",
        "certified/runtests.jl",
        "theory/degree_two.jl",
        "certified/integration.jl",
    )
    for script in scripts:
        log = destination / (script.replace("/", "_") + ".txt")
        print(f"Running {script}; log: {log}", flush=True)
        with log.open("w") as output:
            process = subprocess.Popen(
                ["julia", "--startup-file=no", "--project=.", str(here / script)],
                cwd=root,
                stdout=output,
                stderr=subprocess.STDOUT,
                start_new_session=True,
            )
            try:
                status = process.wait(timeout=180)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
                print(f"Timed out after 180 seconds: {script}", file=sys.stderr)
                return 1
        if status:
            print(f"Failed with status {status}: {script}", file=sys.stderr)
            return 1
        print(f"Passed: {script}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
