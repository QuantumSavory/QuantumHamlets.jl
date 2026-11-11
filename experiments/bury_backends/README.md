# Production BURY backend validation

These measurements exercise the production `backend=:legacy`, `:scan`, and
`:heap` entry points after integration. All retain `Dict{Any,Any}` registry
returns and the original multipartition driver, so their timings differ from
the earlier standalone prototypes with typed registries.

The package suite passes **379,435 assertions**, with **five expected broken
tests** for deferred issues. It compares every backend with a frozen original
implementation, including every labeled simple graph on 1–5 vertices, all
valid targets and explicit seeds on those graphs, multipartitioning, verbose
traces, seed-search output/RNG state, return types, and custom callbacks.
The direct test log is [tests.txt](tests.txt). The standard package test also
passed with the same counts using `Pkg.test(; allow_reresolve=false)` from the
package root; see [pkg_test.txt](pkg_test.txt).

The benchmark script additionally checks ordinary and explicitly seeded
half-target calls on **all 45 GraphML resources**: 90 workload groups with
identical outputs and return types across all three backends, with no graph
mutation. It then measures 12 selected production workloads. No automatic
seed search runs on the large GraphML inputs.

Selected pooled medians, in milliseconds:

| Workload | Legacy | Scan | Heap |
| --- | ---: | ---: | ---: |
| 64×64 grid, ordinary | 9.64 | 10.46 | 0.529 |
| 64×64 grid, explicit seed | 29.35 | 4.57 | 0.635 |
| Initial `barenco_tof_10`, ordinary | 0.204 | 0.215 | 0.0412 |
| Initial `barenco_tof_10`, explicit seed | 0.557 | 0.111 | 0.0532 |
| Evolution `k5`, ordinary | 1.44 | 1.53 | 0.226 |
| Evolution `k5`, explicit seed | 5.48 | 1.05 | 0.359 |
| 32×32 grid, four parts | 1.18 | 1.18 | 0.588 |
| 8×8 grid, automatic seed search | 0.545 | 0.509 | 0.525 |

Heap preserves its advantage on the larger core workloads. The smaller
automatic-search difference is consistent with scoring overhead. Scan remains
useful for small inputs: reduced `barenco_tof_3` with an explicit seed takes
0.959 µs with scan versus 1.542 µs with heap.

Measurements use Julia 1.12.6 on an Apple M1 Pro with one Julia thread.
Three rounds rotate backend order, with `evals=1`, up to 51 samples per trial
(21 for automatic search), and a 0.3-second trial budget. Loading and
compilation are excluded. Each trial has an untimed full GC, and automatic
collection remains enabled. Automatic searches reset the same RNG seed in
untimed setup and include scoring; improvement messages are silenced outside
the timed calls. All 36 saved result rows and 4,964 samples were checked for
consistent counts and medians. These are local measurements, not performance
guarantees for every graph.

From the package root:

```sh
julia --startup-file=no --project=experiments/bury_performance test/runtests.jl
julia --startup-file=no --threads=1 --project=experiments/bury_performance experiments/bury_backends/benchmarks.jl
```

The benchmark command overwrites this directory's results and uses the existing
isolated experiment environment. It imports the earlier measurement helper but
passes production API calls for every competitor. Results are in
[results.tsv](results.tsv), raw timings in [samples.tsv](samples.tsv), and
machine/settings details in [environment.txt](environment.txt).
