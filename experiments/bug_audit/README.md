# Bug audit reproductions

These scripts audit the post-integration code without modifying production
source or the existing package test suite. The [current bug review](../../docs/bug_audit.md)
distinguishes confirmed defects, newly tested scenarios, and documented limits.

From the package root, run the small probes with the existing experiment
environment:

```sh
julia --startup-file=no --project=experiments/bury_performance experiments/bug_audit/known.jl
julia --startup-file=no --project=experiments/bury_performance experiments/bug_audit/api/probes.jl
python3 experiments/bug_audit/costs/run.py
python3 experiments/bug_audit/integration/run.py
```

The numerical sweep passed 127,098 assertions on 7,122 fixtures; API boundary
checks passed 407 assertions; the earlier-findings probes passed 26 observation
assertions. The latter deliberately verify the current erroneous outputs as
evidence. These are standalone audit scripts, not new production regressions
that mark known bugs as acceptable behavior.

| Files | Purpose |
| --- | --- |
| `known.jl`, `known.log` | Current status of score-bound, one-part matching, aliasing, display, and validation findings |
| `api/probes.jl`, `api/probes.log` | Empty seed helpers, failed-call RNG effects, supported input types and conversions |
| `costs/oracles.jl`, `costs/oracles.log` | Independent exhaustive matching and GF(2) rank checks, with cover/witness invariants |
| `costs/README.md` | Numerical coverage and bounded recursion confirmation instructions |
| `costs/recursion_confirm.log` | Fresh-process confirmation of the public matching stack overflow |
| `integration/probes.jl`, `integration/results.log` | Tiny baseline, plotting, display, and script-score probes |
| `integration/disconnected.jl`, `integration/disconnected.log` | Six disconnected baseline cases and the missed zero-cut example |
| `baseline_sha256.json` | Hashes of all 865 pre-existing files before this review |
| `verify_preservation.py` | Verify that those files remain unchanged |

The recursion probe belongs in a disposable subprocess with its supplied
external watchdog, not an interactive research session. The audit stopped the
worker after the confirmed stack overflow. No full resource script or large
automatic seed-search sweep was run.
