# Handoff prompt: QuantumHamlets library and GraphML comparisons

Copy the prompt below into the next agent conversation.

---

Work in `/Users/micciche/Research/QuantumInfo23/JuliaProjects/QuantumHamlets.jl`.
This is a state handoff dated 12 November 2026. Read the current tree before
editing and continue from the user's next instruction. The library changes
and the latest family figures are complete; the latest request was to prepare
this handoff, not to start another implementation pass or experiment.

## User's requests and the immediate unresolved issue

The user initially requested scoped improvements to BURY's exhaustive target
range search and a graph-first near-balanced bipartition API. That work is
implemented and tested, but uncommitted.

The user then requested GraphML comparisons inspired by
`resources/grid_plots.jl` and `resources/qaoa_plots.jl`. The first all-graph
plots were rejected as unattractive and hard to understand. The user steered:

1. Focus on specific circuit families rather than putting all data on a plot.
2. Look for regimes where BURY has an advantage, without biased selection.
3. Compare **BURY, seeded BURY, and METIS**.

Six clearer family figures have now been generated, using three consistent
colors, grouped bars, numerical median labels, and observed min–max whiskers.
All available sizes in the named families are shown. The numerical study
covers every input, including families not selected for the figures, losses,
ties, and invalid outcomes.

The user's two latest clarification questions were:

- **Which seeded BURY, the sqrt one?** Answer: for bipartitions, yes; for
  multipartitions, the existing deterministic minimum-degree residual seeder.
- **Were bipartitions equal, or off by one for odd graphs?** Answer: yes for
  both BURY methods, but **METIS used its default 3% imbalance tolerance**.
  Some METIS partitions therefore differ by more than one vertex. The figures
  flag these cases, but do not enforce identical balance constraints.

**Keep this distinction prominent.** The user has not yet requested a change
to METIS's balance policy or to which seeded algorithm is used. A likely next
discussion is a fairer balance comparison. Do not assume the current METIS
results are strictly balanced, or describe balance-matched subset summaries
as solving the unequal-policy issue. If changing tolerance later, verify actual
returned sizes rather than assuming a tolerance parameter guarantees exact
capacities. Any repair/filter/retry policy should be explicit and produce new
evidence rather than rewriting the existing results.

## Repository state and preservation

Current HEAD is **`a05adb5` (`speedups from codex`)**. No commits were created
in this conversation. There is substantial uncommitted and untracked work.
Run `git status --short`, inspect diffs, and recheck applicable `AGENTS.md`
files. No applicable `AGENTS.md` was found during the work so far.

Preexisting unrelated changes at the start included:

- `.DS_Store`, `docs/.DS_Store` modifications.
- Untracked `docs/bipartition_range_handoff.md` and
  `docs/bipartition_range_review.md`; the review was subsequently updated as
  requested to distinguish implemented work from proposals.
- Untracked `docs/figures/Quantum_Hamlets_JRPG_figures.zip`.

Current library work modifies `src/QuantumHamlets.jl`, `src/api.jl`,
`src/bury.jl`, `test/runtests.jl`, `docs/usage.md`, and `docs/bug_audit.md`,
and adds `test/bipartition_range.jl` and new experiment directories/files.
All GraphML plotting/analysis additions are also untracked. Do not reset,
stash, delete, or overwrite unrelated work. Do not upgrade dependencies.

Preserve these historical/reference assets:

- `test/reference/bury_original.jl` — frozen original, never rewritten.
- `docs/bury_approximation.md`, `experiments/bury_approximation/`.
- Existing benchmark/audit/inventory logs and snapshots, including
  `experiments/graphml_bury/`, `experiments/cost_performance/`, and
  `experiments/bug_audit/`.
- Every existing result directory, including stopped preliminary runs.

The repository `.gitignore` ignores `*.png` and `*.pdf`. Generated PNGs exist
on disk even when they are absent from ordinary git status; SVGs and tables
are available alongside them. No publishing or committing was requested.

## Start with these sources

- `docs/usage.md`, `docs/bug_audit.md`.
- `docs/bipartition_range_review.md` and the earlier
  `docs/bipartition_range_handoff.md` for original requirements.
- `src/QuantumHamlets.jl`, `src/bury.jl`, `src/bury/{heap,scan,legacy}.jl`,
  `src/api.jl`, `src/costs.jl`, `src/partition.jl`, `src/quantum_land.jl`.
- `test/bipartition_range.jl`, `test/bury_validation.jl`,
  `test/bury_backends.jl`, `test/api.jl`.
- **`resources/graphml_family_plots_README.md`** — primary guide to the
  current three-algorithm study, results, figures, caveats, and commands.
- `resources/graphml_family_analysis.jl`,
  `resources/graphml_family_plots.jl`,
  `resources/run_graphml_family_plots.py`, and
  `resources/graphml_family_analysis_validation.jl`.

## Library functionality now implemented

### Public representations and partitioning

`GraphPartition` stores dense membership; constructors and accessors copy
collections. Treat fields as read-only. Graph-first scalar costs include
`matching_cost(g,p)`, `edge_cut_cost(g,p)`, and `pairwise_cut_rank(g,p)`.
Rank is GF(2), and for more than two parts it is a sum of pairwise ranks.

`partition_graph(g,k)` retains its exactly equal-sized multipartition contract:
nonempty graph, integer `k` in `1:n`, and `n % k == 0`. Its seeded algorithm
uses a minimum-degree seed in each residual graph. The multipartition driver
was not rewritten.

New exported helper:

```julia
using Graphs, QuantumHamlets
p = bipartition_graph(path_graph(11); backend=:heap)
p = bipartition_graph(path_graph(11); targets=4:7, objective=pairwise_cut_rank)
```

`bipartition_graph`:

- Uses **ordinary BURY only** and returns `GraphPartition(...; nparts=2)`.
- Defaults to `fld(n,2):cld(n,2)`: one trial for even graphs, both targets for
  odd graphs. Label 2 is the colored side.
- Requires `n >= 2` and a nonempty contiguous ascending integer range in
  `1:n-1`. Rejects floating, reversed, empty, non-unit-step, and empty-sided
  ranges before callbacks.
- Uses graph-first scalar `objective(g,p)`, default matching.
- Evaluates every target in ascending order and keeps the **last** target
  on equal scores. Supports `:heap`, `:scan`, and `:legacy`.
- Calls arbitrary custom objectives once per target; closures are not
  silently inferred to be built-ins.

### Legacy range behavior and optimization

`QuantumHamlets.bury_not_balanced_bipartition(g,min_v,max_v; backend,
bury_v, cost_f)` remains available with its original tuple return and
land-first tuple-cost callback contract.

Preserved behavior:

- Integer ordered bounds in `0:n`; ordinary empty graph with range `0:0`
  remains valid. Seeded empty calls reject before sampling.
- Every target is evaluated, ascending; scores initialize from the first
  candidate, so negative/high/zero scores work; ties retain the last.
- Exact built-in `bury` and `bury_seeded` receive `backend`. A custom
  `bury_v(g,target)` gets exactly two positional arguments and controls its
  own backend. Callback exceptions propagate; no `MethodError` retry.
- Generic paths retain `quantumLand` validation and tuple-cost behavior.
- Returned built-in formats remain `(Dict{Any,Any}, Vector{Bool})`.

Fixed ownership defect: retained generic callback registries and color
vectors are copied. A callback may reuse buffers or mutate them after return
without corrupting the selected result.

Fast path: exact ordinary `bury` plus exact built-in matching/edge-cut/rank
objectives use existing private scalar kernels on Boolean membership.
Heap/scan avoid trial registries, lands, dense integer-label conversion, and
unused edge colors; only the winning return representation is materialized.
Winning membership is copied. Legacy's algorithm body is unchanged, so it
still constructs trial registries internally. Seeded/custom paths keep their
established behavior.

### Algorithm policies not changed

- Heap is the default; scan and legacy remain supported. Valid greedy
  choices, smallest-ID ties, explicit seed behavior, traces, and default-RNG
  consumption are preserved.
- Automatic `bury_seeded(g,target)` calls `best_sqrt_seed`: **ceil(sqrt(n))
  samples with replacement**, default RNG, inner matching objective,
  first-best seed ties, then reruns the selected seed.
- Outer rank/edge-cut range objectives select among matching-selected inner
  candidates. They do not make the inner seed search rank/edge-cut aware.
- `bury_k_seeded` instead uses the deterministic minimum-degree residual
  seeder; it is not the sqrt search.
- Objective-aware seeds, unique/shared seed pools, inner result reuse,
  workspace reuse, and pruning are deferred compatibility decisions.

Do not reintroduce old bounded-score initialization, recursive matching
traversal, or weighted unit-weight matching. Current matching uses native
unweighted cardinality matching and bulk graph construction; rank uses packed
GF(2) elimination. These cost optimizations predated the range implementation.

Do not skip complementary target sizes, infer all targets from one prefix,
assume target-cost monotonicity, or stop at zero. These change results, ties,
callback counts, and/or seeded RNG consumption. Confirmed fixtures include:

- Five vertices, edges `(1,4),(2,3)`: ordinary target 2 colors `{1,5}`, cost 1;
  target 3 colors `{1,4,5}`, cost 0. The near-balanced helper selects sizes 2/3.
- Four vertices, edges `(1,3),(1,4),(2,3),(2,4),(3,4)`: targets 2 and 3 produce
  nonnested sets `{1,2}` and `{1,3,4}`.
- Low-level matching search `0:n` intentionally picks target `n` under last
  ties, giving an empty side. The new helper rejects those endpoints.

## Library validation and range performance

Required full package command passed **459,086 assertions**, including
**1,877** new range tests:

```sh
julia --startup-file=no --project=. -e 'using Pkg; Pkg.test(; allow_reresolve=false)'
```

Log: `experiments/bipartition_range_performance/validation_pkg_test.log`.
Later plotting work did not modify the production library. The new experiment
validation is separate from this full package suite.

`experiments/bipartition_range_performance/README.md` documents 120 bounded
before/after/helper comparisons against the exact old wrapper at `a05adb5`,
using the existing `experiments/bury_performance` environment. All passed.
Built-in wrapper speedups ranged 1.01–5.61x in those fixtures. Generic callback
ownership copies cost 2–12% runtime; six legacy singleton rows increased
allocated bytes, up to 6.6%. Both gains and regressions are retained.
Example: grid225, six-target heap matching, old wrapper 306 us, new wrapper
107 us, helper 92.8 us. These are fixture measurements, not universal promises.

Earlier preserved validation: pre-range full suite 457,209 assertions;
independent cost audit 127,098 assertions on 7,122 fixtures; focused cost
semantics 9,191 assertions. See original documentation for their exact scope.

## GraphML inputs and environments

Collections under `resources/`:

- `graphml_initial`: 20 graphs, 67–652 vertices, 15 odd.
- `graphml_reduced`: 20 graphs, 25–294 vertices, 15 odd.
- `graphml_evolution`: five graphs, 25–1,772 vertices, three odd.

All 45 are connected simple undirected graphs. Node IDs are noncontiguous;
GraphIO maps XML declaration order to Julia vertices `1:n`. Node phase/type/
layout and edge-type attributes are ignored for these unweighted graph-cost
comparisons. Loading was checked against all existing inventory hashes and
independent topology exports under `experiments/graphml_bury/inventory/`.
Do not rerun `inventory.py` casually: it overwrites historical exports.

Run GraphML work from the package root with **`--project=..`**, the existing
parent `JuliaProjects` environment. It has local QuantumHamlets, GraphIO,
EzXML, Metis, CairoMakie, and Graphs installed. This differs from the package
and benchmark environments; do not add their missing research dependencies
to the production package just to run these scripts. No dependencies were
installed or upgraded. Runs used Julia 1.12.6, Metis 1.5.0, GraphIO 0.7.1,
EzXML 1.2.3, CairoMakie 0.15.12, and Graphs 1.14.0 in the parent environment.

The inspiration scripts are not safe imports for a quick check:
`qaoa_plots.jl` launches its sweep at file end, and `sandbox.jl` has top-level
research work. The new modules do not run a sweep when included.

## Current three-algorithm experiment

Primary guide: `resources/graphml_family_plots_README.md`.

Data module: `resources/graphml_family_analysis.jl`.
Plot module: `resources/graphml_family_plots.jl`.
Runner: `resources/run_graphml_family_plots.py`.
Focused tests: `resources/graphml_family_analysis_validation.jl`.

Fixed collection plan, chosen before this expanded run:

- Every one of the 45 graphs, both matching and pairwise GF(2) rank.
- Part counts 2 through 8; every graph at two parts; higher counts only when
  exactly divisible. This yields **79 feasible graph/part-count cases**
  (45 bipartitions + 34 multipartitions) and 236 unsupported combinations.
- Ordinary BURY: deterministic, one result per objective/case. At two parts,
  floor/ceiling targets are selected separately for each objective.
- Seeded BURY, two parts: existing automatic sqrt search through the legacy
  range wrapper, repeated with fixed Julia RNG seeds **0:9**. Reset RNG per
  objective/replicate so the candidate streams agree across objectives.
  Inner selection remains matching even in rank panels.
- Seeded BURY, >=3 parts: existing deterministic minimum-degree residual
  driver, recorded once, not ten duplicated observations.
- METIS: native **KWAY**, edge-cut optimization, equal target weights,
  **UFACTOR=30 (3%)**, fixed native seeds **0:9**. Each membership is reused
  for both score objectives. The wrapper restores caller METIS options.
- Equal integer seeds across Julia and METIS are not paired common random
  choices. Search budgets differ; no runtime or equal-compute claim is made.

Primary figures use medians over nonempty outcomes and observed min–max
ranges, not confidence intervals. `†` marks a case with valid METIS outputs
outside BURY's balance bound; `‡` marks excluded empty-part outcomes.
Raw invalid scores/sizes remain in the tables. An all-invalid comparator is
`unavailable`, not a BURY win. Separate summaries restricted to balance-matched
METIS outputs are a sensitivity check on a selected subset, with counts.

### Authoritative data and final figures

**Frozen data and first draft figures:**
`resources/figs/graphml_families/run_20261112T145558Z/`

- `results.tsv`: **2,706 raw scores**: 158 ordinary BURY, 968 seeded BURY
  (900 sampled bipartition + 68 deterministic multipartition), 1,580 METIS.
- `summary.tsv`: 158 graph/part-count/objective summaries, including medians,
  extrema, total/valid/balanced replicate counts, and both METIS comparators.
- `family_summary.tsv`: 200 complete-family summary rows, including other
  families not selected for figures; lower/tied/higher/unavailable counts.
- `inputs.tsv`, `status.tsv`, `metadata.txt`, `output.log`.
- `validation.log`: **344 focused assertions passed**.
- `audit.txt`: **31,351 independent Python checks passed**, reconstructing
  tables, counts, exclusions, hashes, and seed coverage. This table audit did
  not independently rerun the cost kernels or reconstruct all memberships.

**Final caption-corrected figures, rendered from the same saved data:**
`resources/figs/graphml_families/run_20261112T145758Z/`

Each has PNG and SVG:

1. `multipliers_bipartitions`: initial/reduced x matching/rank, m=4,5,6.
2. `multipliers_multipartitions_matching`: initial/reduced x m=4,5,6,
   all feasible requested counts within 3:8.
3. `multipliers_multipartitions_rank`: same cases, GF(2) rank.
4. `toffoli_bipartitions_matching`: initial/reduced x Toffoli/Barenco,
   all available parameters 3,4,5,10.
5. `toffoli_bipartitions_rank`: same cases, rank.
6. `evolution_bipartitions`: all five stages, both objectives.

All six final PNGs were visually reviewed. SVGs parse, provenance hashes
match, and an independent reviewer checked 117 displayed series. See
`figure_provenance.txt` and `visual_review.txt`. The two seeded policies are
shown in separate bipartition/multipartition figures.

### Findings, including counterevidence

Describe these as fixture-specific **measured cost** patterns, not general
performance guarantees. The multiplier hypothesis was noticed in seed-0 data;
0:9 includes seed 0, so this is exploratory robustness work, not independent
confirmation or a preregistered discovery test.

- **Ordinary BURY on reduced GF multipliers:** all three graphs and all 11
  feasible cases give matching **9 lower / 2 tied / 0 higher**, and rank
  **8 lower / 3 tied / 0 higher**, versus METIS medians.
- **Same-balance examples:** reduced GF(2^4), 132 vertices. At four parts
  (33 each), ordinary matching/rank = 31/31, seeded = 36/34, METIS = 35/35.
  At six parts (22 each), ordinary = 54/51, seeded = 54/50, METIS = 66/63.
  All ten METIS runs meet balance and return those same costs. The other
  nine reduced-GF cases have no balance-matched METIS samples.
- **Counterevidence:** initial GF multipliers give ordinary 2 lower, 1 tied,
  5 higher for each objective. Seeded does not uniformly improve reduced GF:
  at GF(2^5), three parts, minimum-degree seeded costs 65 for both objectives,
  versus ordinary 26 matching / 25 rank and METIS median 30 for both.
- **Sampled seeded, reduced Toffoli/Barenco bipartitions:** across all eight
  graphs, each objective gives **4 lower / 4 tied / 0 higher**. The pattern
  survives the balance-matched subset, but two cases have only 1 or 2 such
  METIS samples (the other six have ten). Multipartition seeded counterexamples
  remain: reduced Barenco parameter 4, three parts, matching 12 vs METIS 8;
  reduced Toffoli parameter 3, three parts, 10 vs 8.
- **Sampled seeded, evolution bipartitions:** matching median lower at all
  five stages; rank 3 lower, 1 tied, 1 higher. Ordinary BURY wins only stage 1.
  Stage 5: ordinary 194 matching / 153 rank; sampled seeded 101/85; METIS
  medians 143.5/118.5 with ranges 124–159 / 101–138. No balance-matched METIS
  sample exists there. At stage 3, seeded matching 27 is above METIS's best
  observed 26 despite beating its median 29; preserve that nuance.

Across all native METIS runs (objectives deduplicated): **779/790 nonempty**,
only **287/779** meeting BURY's balance bound. Empty outcomes are all ten
seeds for evolution stage 1 with five parts, plus seed 3 for reduced GF(2^6)
with seven parts. The latter has nine valid replicates; the former has none.
These produce 22 flagged raw score rows. Do not count invalid results as wins.

### Commands

From the package root:

```sh
# Collect all fixed-seed data and render to a fresh timestamped directory.
python3 resources/run_graphml_family_plots.py

# Prefer this for visual-only revisions: no partitioner rerun.
python3 resources/run_graphml_family_plots.py \
  --input-dir resources/figs/graphml_families/run_20261112T145558Z

# Focused experiment validation.
julia --startup-file=no --project=.. resources/graphml_family_analysis_validation.jl
```

Runner defaults: heap, parts 2:8, seeds 0:9. External limits: 180s startup,
120s per stage, 1,200s overall. It creates a new timestamped output directory,
refuses overwrites, logs all progress, and kills the process group on timeout.
There are no ongoing background runs at this handoff.

## Earlier GraphML outputs retained

`resources/graphml_plots.jl`, `resources/run_graphml_plots.py`,
`resources/graphml_plots_validation.jl`, and
`resources/graphml_plots_README.md` are the earlier two-algorithm workflow.
The README now points to the family guide.

- `resources/figs/graphml_comparison/run_20261112T142507Z/`: preliminary run
  stopped on METIS's empty fifth part in evolution stage 1; preserved.
- `.../run_20261112T142638Z/`: completed two-algorithm, seed-0 data and five
  dense figures, visually superseded by family plots. 45 inputs, 316 raw
  scores, 78 nonempty comparisons and one empty-part case. Focused validation
  passed 1,379 assertions. Do not confuse these results with the later
  2,706-row three-algorithm fixed-seed study.

The older `experiments/graphml_bury/` is a different historical experiment:
backend/prototype performance, not METIS-versus-BURY partition-quality data.

## Next-agent working stance

Respond to the next actual user request. Preserve all current work and data.
Likely follow-ups are stricter balance comparisons, explicit choice of seeded
policy, or further family-specific figures, but none is currently ordered.
Use the saved data for plotting changes. If a policy changes, create a new
clearly labeled experiment and retain the old baseline. Avoid extrapolating
from these fixtures or presenting selected lower-cost cases without their
whole-family losses, ties, capacity constraints, and sampling policy.
