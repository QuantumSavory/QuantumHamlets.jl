# Circuit-family comparisons: BURY, seeded BURY, and METIS

These figures separate named circuit families and graph stages, showing three
algorithms with consistent colors. They replace the crowded all-graph overview
with small panels while retaining the complete numerical results, including
losses, ties, unsupported capacities, and invalid METIS outputs.

The fixed-seed dataset is
[run_20261112T145558Z](figs/graphml_families/run_20261112T145558Z/): **45 graphs**,
**79 feasible graph/part-count cases**, **2,706 raw score rows**, and **158
case/objective summaries**. The full
[results](figs/graphml_families/run_20261112T145558Z/results.tsv),
[case summaries](figs/graphml_families/run_20261112T145558Z/summary.tsv), and
[family summaries](figs/graphml_families/run_20261112T145558Z/family_summary.tsv)
are retained independently of figure revisions.

The final figures are in
[run_20261112T145758Z](figs/graphml_families/run_20261112T145758Z/), rendered
from those saved scores without rerunning the partitioners.

| Figure | PNG | SVG |
| --- | --- | --- |
| Multiplier bipartitions, both objectives | [View](figs/graphml_families/run_20261112T145758Z/multipliers_bipartitions.png) | [Vector](figs/graphml_families/run_20261112T145758Z/multipliers_bipartitions.svg) |
| Multiplier multipartitions, matching | [View](figs/graphml_families/run_20261112T145758Z/multipliers_multipartitions_matching.png) | [Vector](figs/graphml_families/run_20261112T145758Z/multipliers_multipartitions_matching.svg) |
| Multiplier multipartitions, rank | [View](figs/graphml_families/run_20261112T145758Z/multipliers_multipartitions_rank.png) | [Vector](figs/graphml_families/run_20261112T145758Z/multipliers_multipartitions_rank.svg) |
| Toffoli families, matching | [View](figs/graphml_families/run_20261112T145758Z/toffoli_bipartitions_matching.png) | [Vector](figs/graphml_families/run_20261112T145758Z/toffoli_bipartitions_matching.svg) |
| Toffoli families, rank | [View](figs/graphml_families/run_20261112T145758Z/toffoli_bipartitions_rank.png) | [Vector](figs/graphml_families/run_20261112T145758Z/toffoli_bipartitions_rank.svg) |
| Evolution bipartitions, both objectives | [View](figs/graphml_families/run_20261112T145758Z/evolution_bipartitions.png) | [Vector](figs/graphml_families/run_20261112T145758Z/evolution_bipartitions.svg) |

## Observed patterns and counterexamples

**Ordinary BURY is strongest on the supplied reduced GF multipliers.** Across
all three exponents and all 11 feasible part-count cases, its matching cost
is lower than the METIS median in nine cases and tied in two; rank is lower
in eight and tied in three. It is never higher in that subset. Initial GF
graphs show the opposite caution: two lower, one tied, and five higher costs
for both objectives. The pattern concerns these reduced fixtures, not GF
circuits or graph size generally.

Only two reduced-GF cases have any METIS outputs meeting BURY's balance
constraint. All ten native seeds meet it in those cases and return the same
costs. Both preserve ordinary BURY's advantage:

| Reduced GF(2⁴), 132 vertices | Part sizes for all methods | BURY matching / rank | Seeded BURY matching / rank | METIS matching / rank, all ten seeds |
| --- | --- | ---: | ---: | ---: |
| Four parts | `33,33,33,33` | 31 / 31 | 36 / 34 | 35 / 35 |
| Six parts | `22,22,22,22,22,22` | 54 / 51 | 54 / 50 | 66 / 63 |

The other nine reduced-GF cases have **no balance-matched METIS replicate**;
their primary cost differences are comparisons under different capacity
constraints. Seeded BURY does not uniformly improve ordinary BURY: reduced
GF(2⁵) with three parts has seeded cost **65** for both objectives, versus
ordinary matching **26**, ordinary rank **25**, and METIS median **30** for
both. This uses the deterministic minimum-degree seeded variant.

**Sampled seeded BURY improves the evolution bipartitions in this run.** Its
matching median is lower than the METIS median at all five stages; rank is
lower at three, tied at one, and higher at one. Ordinary BURY is lower only
at stage 1 and higher at the other four stages for both objectives. All ten
seeded-BURY results coincide within each stage/objective, an observed
stability over these fixed seeds rather than a guarantee. At stage 3 the
seeded matching score **27** lies inside METIS's observed **26–36** range,
so a lower median does not mean beating every METIS replicate.

Evolution stage 1 provides an equal-capacity example: ordinary and seeded
BURY score matching **3**, rank **2**, while all ten METIS seeds score **5**
for both, with nearest-balanced sizes. At stage 5, seeded BURY scores
matching **101** and rank **85**, versus METIS medians **143.5** and **118.5**;
there are no balance-matched METIS samples there. Balance-matched METIS
samples exist only at stages 1, 2, and 4 (10, 1, and 1 samples respectively).

**Sampled seeded BURY also has a consistent reduced-Toffoli bipartition
pattern.** Across all four reduced Toffoli and all four reduced Barenco
Toffoli graphs, each objective has four lower medians and four ties against
METIS, with no higher costs. Restricting METIS to balance-matched outputs
retains those outcomes. Six cases have ten balance-matched native replicates;
Barenco parameter 10 has two, and Toffoli parameter 10 has one. This result
is specific to the sampled bipartition policy. The deterministic seeded
multipartition variant has counterexamples: reduced Barenco parameter 4 with
three parts scores matching **12** versus METIS median **8**, and reduced
Toffoli parameter 3 with three parts scores **10** versus **8**.

The following counts retain every feasible case within the plotted families,
including multipartitions absent from the bipartition-only Toffoli/evolution
figures. Entries are **lower / tied / higher** than the METIS median over
nonempty runs. They count cases, not independent graphs or statistical wins.

| Family / collection | Feasible cases | BURY matching | BURY rank | Seeded matching | Seeded rank |
| --- | ---: | ---: | ---: | ---: | ---: |
| GF multipliers / initial | 8 | 2 / 1 / 5 | 2 / 1 / 5 | 2 / 2 / 4 | 2 / 2 / 4 |
| GF multipliers / reduced | 11 | 9 / 2 / 0 | 8 / 3 / 0 | 6 / 1 / 4 | 7 / 1 / 3 |
| Toffoli / initial | 6 | 0 / 0 / 6 | 0 / 0 / 6 | 0 / 3 / 3 | 0 / 3 / 3 |
| Toffoli / reduced | 7 | 0 / 1 / 6 | 0 / 2 / 5 | 3 / 2 / 2 | 2 / 3 / 2 |
| Barenco Toffoli / initial | 4 | 0 / 1 / 3 | 0 / 1 / 3 | 2 / 0 / 2 | 2 / 0 / 2 |
| Barenco Toffoli / reduced | 7 | 0 / 0 / 7 | 0 / 1 / 6 | 3 / 2 / 2 | 4 / 1 / 2 |
| Evolution | 8 | 1 / 0 / 6 | 1 / 0 / 6 | 5 / 0 / 2 | 3 / 1 / 3 |

Evolution has one additional `unavailable` comparator: stage 1 with five
parts returns an empty METIS part for every native seed. It is omitted from
the lower/tied/higher totals rather than treated as a favorable BURY result.
All remaining circuit families and the balance-matched sensitivity counts
are in the complete family-summary table linked above.

## Reproduce or replot

Run from the package root using the installed parent `JuliaProjects`
environment:

```sh
# Collect the fixed-seed dataset and render figures into a fresh directory.
python3 resources/run_graphml_family_plots.py

# Render saved results again without rerunning any partitioner.
python3 resources/run_graphml_family_plots.py \
  --input-dir resources/figs/graphml_families/run_20261112T145558Z
```

Defaults are `--backend heap --parts 2,3,4,5,6,7,8` and
`--seeds 0,1,2,3,4,5,6,7,8,9`; both matching and rank are always recorded.
`--output-dir PATH` chooses a fresh destination. The default destination is
a timestamped directory under `resources/figs/graphml_families/`; existing
directories and artifacts are not overwritten. In replot mode, the saved
rows determine the data, seeds, and algorithms: collection options do not
change those results.

The runner invokes Julia with `--startup-file=no --threads=1 --project=..`
without installing or changing dependencies. Its external limits are 180
seconds for loading, 120 seconds per stage, and 1,200 seconds overall. Full
progress is captured in `output.log`, with throttled terminal output. Failures
preserve partial evidence. Including either Julia source defines its module
without launching collection or rendering.

## Comparison protocol

All `.graphml` files directly within `resources/graphml*` are included: 20
initial, 20 reduced, and five evolution graphs. Filenames determine family
and circuit parameter. The topology is simple, undirected, and unweighted;
node phase/type/layout and edge-type attributes do not enter the objectives.
Input hashes and topology sizes are in
[inputs.tsv](figs/graphml_families/run_20261112T145558Z/inputs.tsv).

Every graph supports a near-balanced bipartition. For three through eight
parts, a case is feasible only when the vertex count is divisible by the
part count. There are 45 bipartition cases, 34 feasible multipartition cases,
and 236 unsupported graph/part-count combinations. Unsupported cases remain
in [status.tsv](figs/graphml_families/run_20261112T145558Z/status.tsv); they
are not padded or assigned a different BURY capacity policy.

| Algorithm | Two parts | Three or more parts | Observations per objective/case |
| --- | --- | --- | --- |
| BURY | Ordinary `bipartition_graph`, selecting between floor/ceiling targets separately for each objective | Ordinary equal-sized `partition_graph` | One deterministic result |
| Seeded BURY | Existing sampled seed search at each floor/ceiling target, followed by outer objective selection | Existing deterministic minimum-degree seed on each residual graph | Ten bipartition replicates; one deterministic multipartition result |
| METIS | Native KWAY, equal target weights, edge-cut objective, `UFACTOR=30` | Same policy | Ten native seed replicates |

The ten fixed seed values **0 through 9** were chosen before this expanded
run. For sampled BURY, `Random.seed!(seed)` is reset separately for each
objective and replicate. Its inner sampled search **always selects by
matching cost**, including in rank panels; the outer target search then
minimizes the selected objective among those matching-selected candidates.
The established sampled-seed order, replacement sampling, first-best inner
ties, and last-best target ties remain unchanged. The minimum-degree
multipartition variant does not sample seeds and must not be interpreted as
ten stochastic trials.

METIS runs once per graph, part count, and native seed; each resulting
membership is scored by both objectives. Its KWAY search optimizes edge cut,
not matching or rank. The same integer seed in Julia and native METIS does
**not** provide paired random choices between algorithms. These are cost
comparisons under the documented search policies, with no equal-compute-budget
or runtime-efficiency claim.

All algorithms are scored by the same scalar functions:

- `matching_cost`: sum of maximum matching sizes across unordered part pairs.
- `pairwise_cut_rank`: sum of GF(2) biadjacency ranks across unordered part
  pairs. For more than two parts, this is a sum of pairwise ranks, not one
  bipartite cut rank.

Lower scores are better. Ordinary BURY and seeded multipartition results
are recorded once; duplicate deterministic observations are not introduced.
The 2,706 rows comprise **158 ordinary BURY**, **968 seeded BURY** (900 sampled
bipartition and 68 deterministic multipartition), and **1,580 METIS** scores.

## Balance, failures, and uncertainty

The primary figures summarize every **nonempty** METIS output under its
standard 3% KWAY tolerance. BURY requires exact part sizes for multipartitions
and floor/ceiling sizes for bipartitions. These capacity constraints can differ.
Actual sizes and the `matches_bury_balance` flag are recorded for every result.

An empty-part METIS output is retained in the raw table with `nonempty=false`
and flagged in `status.tsv`. It contributes neither a plotted score nor a
summary median/minimum/maximum. The number of valid replicates remains
explicit. An all-invalid METIS case has no comparator; it is not counted as
a BURY advantage. The valid BURY results remain available.

In this dataset, **779 of 790 native METIS outcomes are nonempty**, and **287
of those 779** meet BURY's balance requirement. The 11 empty outcomes are all
ten seeds for evolution stage 1 with five parts, plus seed 3 for reduced
GF(2⁶) with seven parts. Each outcome is scored twice, yielding **22 flagged
raw score rows**. Of 315 requested graph/part-count status rows, 236 are
unsupported, 77 have every METIS replicate nonempty, and two have at least
one empty-part result. The partially valid reduced-GF case contributes nine
valid METIS scores per objective to its summaries.

`summary.tsv` provides two METIS comparisons for every case:

1. All nonempty native-seed outcomes, used in the main figures.
2. Only outcomes meeting BURY's balance constraint, reported in
   `metis_balanced_median/min/max` and the balance-matched comparison columns.

The second is a sensitivity check on a potentially smaller selected subset,
not a separate independent experiment. A case with no balance-matched outcome
is marked `unavailable`. Both total/valid replicate counts and exact-balance
counts are retained; `family_summary.tsv` counts lower, tied, higher, and
unavailable cases for every family under both definitions.

Bars show medians and whiskers show **observed min–max ranges across up to
ten valid replicates**. Deterministic variants have one result and no replicate
range. These are not confidence intervals, significance tests, or estimates
from independent graph samples. Graph/part-count cases share source graphs.
Medians include all valid fixed seeds rather than choosing the best seed.

## Family figures

Six figures are rendered in PNG and SVG:

| File stem | Panels and coverage |
| --- | --- |
| `multipliers_bipartitions` | Initial/reduced rows × matching/rank columns; every GF(2ᵐ) multiplier exponent, `m=4,5,6`, at two parts |
| `multipliers_multipartitions_matching` | Initial/reduced rows × exponents 4/5/6; every feasible part count from 3 through 8 |
| `multipliers_multipartitions_rank` | The same complete multipartition cases, scored by rank |
| `toffoli_bipartitions_matching` | Initial/reduced rows × Toffoli/Barenco families; all available circuit parameters at two parts |
| `toffoli_bipartitions_rank` | The same complete Toffoli and Barenco bipartition cases, scored by rank |
| `evolution_bipartitions` | All five evolution stages, with matching and rank side by side |

Blue indicates ordinary BURY, green seeded BURY, and orange METIS. A dagger
(`†`) flags a case where some valid METIS runs have part sizes outside BURY's
constraint. A double dagger (`‡`) flags excluded empty-part runs. Unsupported
part counts are omitted and values are not interpolated. All other circuit families and
every feasible multipartition case remain in the data and summary tables,
including cases outside these selected figure panels.

The multiplier pattern was first noticed in the earlier seed-0 comparison.
This ten-seed analysis is therefore an **exploratory robustness follow-up**,
not independent confirmation or a preregistered discovery test. Retaining
initial multiplier graphs, complete family sizes, Toffoli/Barenco comparisons,
and evolution stages supplies counterevidence alongside favorable examples.
Any advantage should be described as an observed difference on these supplied
fixtures, not a general graph-size, graph-density, or circuit-family guarantee.

## Evidence and validation

| File | Purpose |
| --- | --- |
| `results.tsv` | All raw scores, sizes, algorithm variants, replicate IDs, seeds, and validity/balance flags |
| `summary.tsv` | All case/objective summaries, including valid counts and both METIS comparator definitions |
| `family_summary.tsv` | Complete family counts, retaining lower/tied/higher/unavailable outcomes |
| `inputs.tsv`, `status.tsv` | Input hashes/topology sizes and requested-case coverage/failures |
| `metadata.txt` | Algorithm policies, package versions, implementation hashes, and collection elapsed time |
| `figure_provenance.txt` | The input results path/hash and plotting source hash for each rendering |
| `output.log` | Exact command, watchdog limits, full progress, and completion status |

```sh
julia --startup-file=no --project=.. resources/graphml_family_analysis_validation.jl
```

Focused validation passed **344 assertions in 10.2 seconds** using tiny
fixtures; see [validation.log](figs/graphml_families/run_20261112T145558Z/validation.log).
It checks seeded membership and subsequent RNG draws against an
explicit target loop, deterministic multipartition variants, objective
scoring, invalid-input and overwrite preservation, and summaries containing
empty, off-balance, and all-invalid METIS outcomes. The data module does not
change production partitioners or scoring kernels.

An [independent completed-data audit](figs/graphml_families/run_20261112T145558Z/audit.txt)
passed **31,351 checks**, covering source hashes for all 45 graphs, raw-row
coverage, summary reconstruction, invalid outcomes, and family counts. The
final [figure provenance](figs/graphml_families/run_20261112T145758Z/figure_provenance.txt)
records the saved-data hash and plotting-source hash. The initial data/figure
run and later caption-corrected rendering are both preserved.
