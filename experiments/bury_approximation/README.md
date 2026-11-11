# BURY approximation proofs and experiments

These additions investigate maximum matching cost for an exactly balanced bipartition. The [mathematical report](../../docs/bury_approximation.md) gives the full proofs and interpretation.

The principal result is an infinite connected planar family of maximum degree four with optimum one and BURY matching cost at least `ceil(n/6)`, regardless of seed or vertex labeling. The separator prototype gives rigorous bounds on individual instances and repairs all nonempty optimum-zero/one inputs with a search cap of one. Use cap zero for the empty graph. It does not give a general polynomial-time constant approximation. A simpler component guard gives a factor two on graphs of maximum degree at most two.

## Reproduce the main checks

Run from the package root with the existing Julia environment:

```sh
julia --startup-file=no --project=. experiments/bury_approximation/counterexamples/path_hubs.jl
julia --startup-file=no --project=. experiments/bury_approximation/certified/runtests.jl
julia --startup-file=no --project=. experiments/bury_approximation/theory/degree_two.jl
julia --startup-file=no --project=. experiments/bury_approximation/certified/integration.jl
```

Alternatively, `python3 experiments/bury_approximation/run_checks.py` runs those four commands sequentially, with a 180-second external limit per process, and writes fresh output to `rerun_logs/`. It stops at the first failure or timeout and does not overwrite the original evidence logs.

All inputs are bounded synthetic cases. None of these commands runs the GraphML datasets or full research sweeps. The original runs used Julia 1.12.6 and passed **374,323 assertions** across the four main checks:

| Script | Saved result | Assertions |
| --- | --- | ---: |
| `counterexamples/path_hubs.jl` | `counterexamples/path_hubs_results.txt` | 12,861 |
| `certified/runtests.jl` | `certified/results.txt` | 158,654 |
| `theory/degree_two.jl` | `theory/degree_two_results.txt` | 202,652 |
| `certified/integration.jl` | `certified/integration.txt` | 156 |

The theorem is an analytical statement about all family sizes; these finite tests confirm the relevant behavior of the current implementation.

## Contents

- `counterexamples/path_hubs.jl`: the strongest construction, with `path_hub_obstruction(r; b=4)` for `r≥3`. It tests every seed on each instance and on two random relabelings, plus backend equivalence checks.
- `certified/certified_bury.jl`: independent module exposing `certified_bury` and `component_guarded_bury`. Both call the existing production BURY; neither changes its greedy rule. The first returns a partition, exact matching cost, lower bound, displayed factor, and optimality status. The second returns a partition after checking whether whole components can give cost zero.
- `certified/runtests.jl`: exhaustive tiny-instance oracle based on direct enumeration of cuts and matchings, independently of the separator characterization and the package's matching implementation.
- `certified/integration.jl`: zero-cost packing regressions, validation, fallback equivalence, and optimal repair of the main counterexamples.
- `certified/proof.txt`: proof and complexity of the separator characterization and the component guard.
- `theory/degree_two.jl`: exhaustive checks of BURY's at-most-two-crossing-edges invariant for graphs of maximum degree two, including all targets and explicit seeds.
- `figure.jl`: generates the 30-vertex counterexample illustration from the production partition and an explicit optimum. Outputs are `figures/bury_obstruction.png` and `.pdf`.

Earlier exploratory constructions are retained in `counterexamples/bounded_degree.jl` and `counterexamples/all_seeds.jl`, with their corresponding saved logs. The path-hub construction supersedes them mathematically. `theory/separator_oracle.py` is an additional independent small-graph check of the separator equivalence; its counts are not included in the main total above.

## Scope of changes

This study adds only this experiment directory and `docs/bury_approximation.md`. It makes no edits to existing source files, package tests, dependencies, or historical benchmark results. The certificate cap must be chosen explicitly for large workloads: enumeration is polynomial only for a fixed cap, with a potentially high exponent. No universal constant-factor claim is made for a fixed cap on general graphs.
