# BURY performance experiment

This directory preserves the standalone prototypes and measurements from the
study before production backend integration. The package now provides
`backend=:heap` (default), `:scan`, and `:legacy`; these experimental modules
remain independent. Benchmark rows labeled `current` mean the original
implementation, and the scripts now request `backend=:legacy` explicitly so
reruns do not silently use the new heap default.

Recorded timings and reports describe the original experiment, including its
typed prototype registries. Production wrappers preserve the original public
return types, so use separate production measurements when assessing those
entry points.

From the package root, prepare the isolated benchmark environment:

```sh
julia --startup-file=no --project=experiments/bury_performance -e 'using Pkg; Pkg.instantiate()'
```

Its manifest retains the dependency versions from the package's existing manifest.
The local `QuantumHamlets` dependency points to `../..`. BenchmarkTools is a direct
dependency of this experiment only.

Run exact-output comparisons before measuring:

```sh
julia --startup-file=no --project=experiments/bury_performance experiments/bury_performance/runtests.jl
```

Then run benchmarks without other tests or benchmarks executing concurrently:

```sh
julia --startup-file=no --project=experiments/bury_performance experiments/bury_performance/benchmarks.jl
```

After that process finishes, collect additional samples for the largest cases:

```sh
julia --startup-file=no --project=experiments/bury_performance experiments/bury_performance/slow_cases.jl
```

This uses four seconds per trial and writes `slow_results.tsv` and
`slow_samples.tsv`. The report uses these better-sampled results for the 128×128
grid, while retaining the original broad-run results for inspection.

This overwrites `results.tsv`, `samples.tsv`, `environment.txt`, and `profiles.txt`
inside this directory. Results compare warmed calls, include output construction,
and exclude graph generation and loading. Three rounds rotate implementation order.
All calls use one evaluation per sample, at most 101 samples and 0.4 seconds per
trial, with the default RNG reset outside timing. Slow calls have fewer samples;
the actual counts and quartiles are saved. Allocation figures are per-call
allocation traffic, not peak memory.

| File | Role |
| --- | --- |
| `scan.jl` | Direct minimum scans, incremental counts, and fewer temporaries |
| `heap.jl` | Indexed heaps preserving minimum-weight and vertex-ID tie order |
| `multi.jl` | Color-only multipartition driver and direct final-part assignment |
| `search.jl` | Existing seed/range search policy with experimental kernels and typed registries |
| `runtests.jl` | Exact registry/color comparisons, boundary targets, and graph preservation |
| `benchmarks.jl` | Synthetic and bundled QAOA workloads, raw timings, allocations, and profiles |
| `results.tsv` | Median, quartiles, minimum, allocations, and actual sample counts |
| `samples.tsv` | Individual timing samples and round IDs |
| `slow_cases.jl` | Longer trials for the largest ordinary and seeded grid cases |
| `baseline_sha256.json` | Hashes of all files present before this experiment |

The prototypes intentionally return typed `Dict{Int,Int}` registries. Equivalence
means identical registry entries and Boolean color vectors on valid inputs, not
identical dictionary element types/iteration order, diagnostic logging, or
invalid-input behavior. No production bug fixes are included. The search adapter
retains the old matching evaluator and candidate policies, including repeated
sampled seeds, score ties, and the range search's existing initial bound.

The preservation snapshot records the files before this historical experiment:

```sh
python3 experiments/bury_performance/verify_preservation.py
```

It is expected to report differences after the authorized production backend
integration. The historical hashes and recorded results are retained unchanged;
the script is not a current package test. Run `Pkg.test()` for production backend
compatibility checks.
