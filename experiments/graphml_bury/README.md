# GraphML BURY backend experiment

This historical experiment compares the original BURY implementation (`legacy`)
with the unchanged `scan` and `heap` prototypes from `../bury_performance/`. It covers the
45 files in `resources/graphml_initial`, `graphml_reduced`, and
`graphml_evolution`. The package now provides production backend options; the
worker explicitly uses `backend=:legacy` for its reference so the heap default
does not change the meaning of a rerun. Recorded data still describe the
prototype entry points, not the production wrappers.

See the [findings report](../../docs/graphml_backend_findings.md) and
[per-graph guide](per_graph_guide.md). Raw properties and timings are retained
so that the recommendations can be checked against individual files.

## Reproduce

From the package root, use the existing isolated experiment environment:

```sh
julia --startup-file=no --project=experiments/bury_performance -e 'using Pkg; Pkg.instantiate()'
python3 experiments/graphml_bury/inventory/inventory.py
python3 experiments/graphml_bury/run_benchmarks.py > experiments/graphml_bury/run.log 2>&1
python3 experiments/graphml_bury/analyze.py > experiments/graphml_bury/analysis.txt # Historical preservation check; see below
python3 experiments/graphml_bury/make_guide.py
julia --startup-file=no --project=experiments/bury_performance experiments/graphml_bury/plot_results.jl
```

The commands overwrite generated results in this experiment directory. Run
measurements serially, without other tests or benchmarks. The prototype files,
package source, and package environment are read only. Plot generation is
separate from measurement.

`analyze.py` recomputes comparisons and validates the raw measurements, then
checks the pre-integration file hashes. That final preservation check is
expected to fail after the authorized production changes. The historical
snapshot is intentionally unchanged; use `Pkg.test()` for current backend
compatibility. The saved `analysis.txt` records the successful original audit.

## Workloads and bounds

For each graph, the primary target is `floor(n/2)`. Both ordinary BURY and
explicitly seeded BURY are measured; the explicit seed is vertex 1 after mapping
GraphML node declaration order to `1:n`. Automatic seed selection is a separate
workflow with a fixed RNG seed, the original sampled-seed policy, and the same
matching evaluator. The `seed` argument stored in the plan/results is ignored
for ordinary and automatic workloads.

Additional quarter and three-quarter targets are planned for initial/reduced
`barenco_tof_3`, `barenco_tof_10`, `gf2^6_mult`, and `qcla_com_7`, plus evolution
`k1`, `k3`, and `k5`. These are single bipartitions; no multipartition driver or
exhaustive seed search is benchmarked. Odd graph sizes leave the two halves
differing by one vertex.

The controller starts each circuit family at its smallest parameter and applies
the following bounds:

- A warmed preflight call over **one second** stops the rest of that graph and
  prevents escalation to larger parameters of that family.
- An external **eight-second stage watchdog** kills the Julia process group if
  a preflight, trial, or other job stage fails to return. It also covers native
  matching calls. The Julia launcher and its child are in the same new process
  group, so both are stopped.
- Worker loading/compilation has a separate **120-second startup allowance**.
- A **600-second soft run budget** is checked between jobs, after initial
  startup. Active stages retain their eight-second hard limit.

The final run completed 174 groups, stopped work on `k5` after its legacy automatic
preflight exceeded one second, and skipped its four subsequent target variants.
No stage timeout fired. The half-target ordinary and explicitly seeded `k5`
comparisons had already finished. Its skipped core calls were not demonstrated
to be slow; the graph-level stop rule is deliberately conservative.

## Timing and correctness

All competitors are compiled on a small graph first, then checked on the real
input before timing. Exact registry entries and color flags must agree with
legacy; graph nonmutation, colored count, and registry consistency are checked.
Automatic searches must also leave the next RNG draw identical.

Three rounds rotate backend order. BenchmarkTools uses `evals=1`, a maximum of
101 samples per trial (11 for automatic search), and a 0.2-second trial budget.
Medians and quartiles pool the actual samples across rounds. A full GC occurs
once per workload group outside timing; `gctrial=false` avoids nine full
collections for tiny workloads, and automatic garbage collection stays enabled.
Graph loading and compilation are excluded; output construction is included.
Automatic timing includes seed search, matching scoring, and the final BURY
call. All search improvement messages are silenced outside timed calls.

Trials can have few samples when one evaluation exceeds the trial budget.
Evolution `k4` automatic has only three samples per backend, so its small
heap/scan difference should not be ranked confidently. The slow `k5` automatic
entry is a single preflight observation, not a benchmark median; scan and heap
were not attempted for that workflow.

Both prototypes return `Dict{Int,Int}`, while legacy returns `Dict{Any,Any}`.
These are comparisons of the experimental entry points, including
that output-construction difference. They establish partition-content
equivalence on tested inputs, not identical logging, dictionary iteration,
invalid-input handling, or complete API equivalence.

## Artifacts

| File | Contents |
| --- | --- |
| `inventory/inventory.tsv` | All graph properties, hashes, families, topology/map paths |
| `inventory/README.md` | Property-based expectations recorded before timing |
| `inventory/inventory.py` | XML parsing, validation, degree/component statistics, topology export |
| `plan.json`, `status.tsv` | Planned jobs and exact completed/stopped/skipped status |
| `worker.jl`, `run_benchmarks.py` | Julia measurements and external watchdog/controller |
| `results.tsv`, `samples.tsv` | 522 backend result rows and 121,851 raw samples |
| `comparison.tsv` | Per-workload times, ratios, and backends within 10% of the lowest median |
| `validation.tsv` | Exact-output, graph, balance, and RNG checks |
| `jobs/` | Per-workload source records, including the stopped preflight |
| `analyze.py`, `analysis.txt` | Recomputed medians/counts and preservation checks |
| `backend_comparison.png`, `.pdf` | Scan/heap comparison for all 45 primary graph cases |
| `environment.txt` | Julia, dependency, CPU, thread, and benchmark settings |
| `baseline_sha256.json` | Hashes of all 210 files present before this experiment |

The inventory parses XML attributes rather than relying on the earlier QAOA
regex reader. Node/edge circuit metadata is recorded but excluded from the
simple unweighted topology used by BURY. Source GraphML IDs are noncontiguous;
the exported node maps make the relabeling explicit.
