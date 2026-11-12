# QuantumHamlets versus min-cutrank

**Results summary — 12 November 2026**  
Saved experiment: `run_20261112T161827Z`

Sampled sqrt-seeded BURY produced lower median GF(2) cut rank than random-start
min-cutrank on **32 of 43 circuit graphs with complete native annealing results**,
with six ties and five losses. Initializing the annealer with sqrt BURY and
retaining its best observed partition gave lower median rank on **40 of all 45
circuit graphs by 30 seconds**, with four ties and one loss. Ordinary BURY alone
was less consistent: its native rank was lower on 15 graphs and higher on 25.

These findings concern exactly balanced bipartitions and a **tested ten-sweep
annealing schedule**. A separate 100-sweep sensitivity study changed some method
rankings. The evidence supports sqrt BURY as a useful standalone heuristic and
annealing initializer on this corpus; it does not establish superiority over
fully tuned min-cutrank or optimality of the returned cuts.

The [full illustrated report (PDF)][pdf] and [editable report (DOCX)][docx] contain
additional tables, figures, and the complete circuit appendix.

## Comparison design

The main study includes **67 graphs**: all 45 supplied GraphML circuit graphs
(20 initial, 20 reduced, five evolution stages), four square grids with sides
4, 6, 8, and 10, and 18 Erdős–Rényi graphs. The latter use vertex counts
40, 80, and 120; edge probability `p=c/n` for `c=2,5`; and graph seeds
1000–1002 at each size and density. Disconnected samples are retained. Reusing
those seeds creates related fixtures: at a fixed size and seed, the `c=2` graph
is a subgraph of the `c=5` graph. These are not 18 independent population samples.
GraphML attributes such as phase, type, and layout are ignored; every method
receives the same saved unweighted topology.

**Every partition has capacities `floor(n/2)` and `ceil(n/2)`, verified from its
actual membership.** The primary objective is the GF(2) rank of the adjacency
matrix across the cut, equal to QuantumHamlets' `pairwise_cut_rank` for two parts.
This study does not adapt the annealer to multipartitioning.

| Method | Construction |
| --- | --- |
| Ordinary BURY | Heap BURY; evaluate floor/ceiling target candidates and select by rank. |
| Sqrt BURY | Existing automatic seeder: `ceil(sqrt(n))` samples with replacement, inner seed selection by **matching cost**, outer target selection by rank. |
| Random-start SA | Upstream min-cutrank annealer, initialized with a random fixed-size partition. |
| BURY + SA | Same annealer, initialized with the ordinary BURY partition. |
| Sqrt BURY + SA | Same annealer, initialized with that replicate's sqrt BURY partition. |

The rank objective selects BURY candidates; it does not change the inner greedy
construction. Sqrt seeding retains first-best inner seed ties and reruns the
winning seed. It is distinct from the minimum-degree residual seeder used by
the existing multipartition driver. On odd graphs, BURY tries both target sizes
and keeps the last on equal outer scores. Before annealing, membership is
complemented if necessary to put exactly `floor(n/2)` vertices on the row side;
this preserves all cut costs and fixes the annealer's orientation policy.

The unchanged upstream `cut_rank_annealing_row_formula` was pinned to commit
[`86c1d8f`][upstream]. The main schedule uses ten linearly spaced temperatures
from 1.0 to 0.1, following its example/test schedule. This is an explicit schedule
choice, not a universal library default. Algorithm seeds are 0–9; random initial
partitions use separate Python seeds 100000–100009. Equal integer seeds across
Julia and Python do not represent paired common random choices. Ordinary BURY
has one deterministic result per graph; stochastic methods have ten replicates.

Each solve has a 30-second elapsed-time cap and a 75-second external watchdog.
There is one invocation per replicate, without restarts to fill unused time.
Two output views are kept separate:

- **Native:** the final state after all ten sweeps. A comparison is unavailable
  unless all planned replicates complete. Timeouts are never counted as losses.
- **Incumbent by 30 seconds:** the best rank observed by that timestamp, including
  the initializer. A transparent observer retains completed candidate states
  without changing proposals or RNG draws. This is an explicit wrapper policy;
  it is not the upstream final-state return policy.

There are 737 BURY baseline results and 2,010 annealing runs: **1,950 completed,
60 capped, and zero execution errors**. All caps occur at evolution stages 4
and 5, across all three annealing initializers and ten seeds. Valid incumbents
are available for every method at 30 seconds. Of the completed annealing runs,
**236 finish with worse rank than an earlier incumbent**, making the distinction
between these output views consequential.

## Main circuit results

Entries count graphs with **lower / tied / higher / unavailable median rank**
relative to random-start SA in the same output view. Lower is better; ordinary
BURY's single deterministic value is compared with the reference median.

| Method | Native output, 45 graphs | Incumbent by 30 seconds, 45 graphs |
| --- | --- | --- |
| Ordinary BURY | 15 / 3 / 25 / 2 | 16 / 3 / 26 / 0 |
| Sqrt BURY | 32 / 6 / 5 / 2 | 34 / 6 / 5 / 0 |
| BURY + SA | 26 / 10 / 7 / 2 | 29 / 8 / 8 / 0 |
| Sqrt BURY + SA | 38 / 3 / 2 / 2 | 40 / 4 / 1 / 0 |

The initial and reduced collections behave differently. Sqrt BURY has lower
native median rank on **all 20 initial circuits**. On the 20 reduced circuits,
it has **11 lower, five tied, and four higher** medians. Ordinary BURY has
13 lower / one tied / six higher on initial circuits, but only two lower /
one tied / 17 higher on reduced circuits.

Whole-family results retain counterexamples. Across the six initial/reduced GF
multipliers, sqrt BURY has five lower native medians and one tie. Across eight
Toffoli graphs it has six lower and two higher; across eight Barenco graphs,
five lower and three tied. Ordinary BURY is higher on all eight Toffoli graphs.
The [remaining-family figure][remaining] covers all other circuits rather than
discarding families with less favorable outcomes.

![Native cut rank for all multiplier, Toffoli, and Barenco parameters][native]

*Bars show run medians; whiskers show observed minimum–maximum ranges, not
confidence intervals. All available parameters in these families are included.*

Standalone sqrt BURY's five losses against the 30-second random-start incumbent
are reduced `rc_adder_6` (6 versus 3.5), `tof_3` (4.5 versus 2), `tof_4`
(4 versus 3.5), `vbe_adder_3` (5 versus 3), and evolution stage 2
(12 versus 11.5). The sole sqrt-hybrid loss is reduced `vbe_adder_3`
(3.5 versus 3).

## BURY as an annealing initializer

The sqrt-initialized incumbent improves on its own initializer in **186 of 450
circuit replicates**, tying it in 264. At the graph-median level, 21 of 45 improve
and 24 tie. Retaining the initializer prevents regression by construction;
the strict improvements are the substantive evidence for added value. Native
final states can regress: 71 of 430 completed circuit sqrt-hybrid runs finish
above their initializer, by up to four rank units.

The evolution series illustrates both the benefit and the limits. These are
median **incumbent ranks available by 30 seconds**:

| Stage | Vertices | BURY | Sqrt BURY | Random SA | BURY + SA | Sqrt BURY + SA |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 25 | 2 | 2 | 2 | 2 | 2 |
| 2 | 159 | 14 | 12 | 11.5 | 11 | 6 |
| 3 | 562 | 46 | 25 | 33 | 29 | 20 |
| 4 | 1,159 | 99 | 65 | 96.5 | 64 | 39 |
| 5 | 1,772 | 153 | 85 | 204 | 128.5 | 59.5 |

No native annealing result completes at stages 4 or 5. Their incumbent results
measure bounded search quality, not convergence. See the [evolution figure][evolution]
and [hybrid family figure][hybrid].

## Synthetic controls and schedule sensitivity

All methods' native outputs on the four grids have ranks 4, 6, 8, and 10,
respectively; larger-grid optimality is not assumed. Across the 18 sparse random
graphs, sqrt BURY has seven lower / eight tied / three higher native medians
than random SA. The sqrt hybrid has 12 lower / six tied / zero higher under
both native and 30-second views. The [synthetic figure][synthetic] shows all
graph seeds. Its ER bars are medians of three graphwise run medians, with
whiskers spanning those graphwise medians; they are not confidence intervals.

A separately predeclared **100-sweep supplement** uses 11 structural
representatives, seeds 0–2, the same three initializers, and the same cap.
Of 99 runs, 81 complete and 18 are capped, with no errors. All capped cases
are initial GF(2^6) and evolution stage 5. Comparing identical graph,
initializer, and seed keys against the ten-sweep study gives:

| Initializer | 100-sweep incumbent lower / tied / higher at 30 seconds (33 pairs) |
| --- | --- |
| Random | 18 / 13 / 2 |
| Ordinary BURY | 9 / 22 / 2 |
| Sqrt BURY | 6 / 25 / 2 |

Changing to 100 temperatures slows cooling and changes the trajectory; this is
not a continuation of the ten-sweep run. On the selected 11 graphs, the
sqrt-hybrid versus random-SA graph-median comparison changes from
**seven lower / four tied / zero higher** to **six lower / three tied / two
higher**. Both new losses are `n=120, c=2` graphs, seeds 1001 and 1002.
This schedule dependence limits any general claim from the main results.
The [sensitivity summary][sensitivity] retains all budgets and native exclusions.

## Secondary costs, timing, and interpretation

Matching and edge cut are measured on the same rank-selected memberships;
they are not separate objective-specific optimizations. Across 45 circuits at
30 seconds, sqrt BURY has lower / tied / higher matching medians on 39 / 4 / 2
graphs, and edge-cut medians on 40 / 1 / 4. The sqrt hybrid gives 42 / 3 / 0 for
matching and 40 / 1 / 4 for edge cut. Tradeoffs remain: on initial `tof_3`,
the sqrt hybrid has rank 3 versus random SA's 4, but edge cut 12 versus 7.5.

Separate three-repeat warmed Julia measurements give medians across
graph/method/seed timing medians of approximately **10.6 microseconds for
ordinary BURY** and **0.616 milliseconds for sqrt BURY**. These were collected
on an Apple M1 Pro with 16 GiB RAM, macOS 14.8.9, and one Julia thread. Main
annealing collection allowed two Python workers. Hybrid budget accounting adds
the original measured Julia initializer time to the measured Python phase;
it is a composed cost, not a fresh integrated pipeline measurement. Python
timing includes dense setup and observation overhead, while imports, input
file reads, and post-solve validation are excluded. Runtime, language, and
contention differences preclude a precise language-independent speedup claim.

![Incumbent rank against elapsed solve time on predetermined representatives][time]

*Curves begin when all required runs have candidates. Completed solves retain
their incumbent at later budgets; no extra searching or restarts are implied.*

The experiment independently checks actual capacities and GF(2) ranks, with
all **3,563 distinct published main memberships** rescored in both Python and
Julia for rank, matching, and edge cut. Upstream swap updates, observer
transparency, time-budget eligibility, and tiny exact optima were also audited.
Tiny-graph optima do not certify the large instances. Related circuit sizes,
initial/reduced pairs, and evolution stages are not independent population
samples; no population significance or universal performance claim is made.

**METIS is outside this strict-capacity comparison.** The earlier GraphML study
used its default 3% imbalance tolerance. Filtering that earlier dataset to
balanced outputs would not establish an equal-policy comparison.

## Data, figures, and reproduction

The saved [protocol][protocol] records the input plan, seeds, schedules, budgets,
and provenance. The [complete study archive][archive] contains the experiment
README with commands, environments, validation details, and all six final
figure designs in PNG/SVG/PDF. During this study, production library code was
unchanged and no dependencies were installed or upgraded.

- [Per-graph summaries][summary] and [complete family summaries][families].
- [Input manifest][inputs] and [schedule sensitivity summary][sensitivity].
- The archive preserves `results.tsv`, `memberships.tsv`, individual run records,
  independent comparison audits, table/rescoring validation, and reproduction
  code. The 13,735 score rows include multiple views and budgets of shared runs;
  they are not 13,735 independent trials.

The final reports, six PNG/SVG figure pairs, and selected summary tables remain
directly accessible under `docs/results/min_cutrank/`. In the archive, final
figures retain their original `figures_v2/` paths; the earlier `figures/`
directory and other drafts are also preserved. Follow the archive index's
restoration instructions before using the historical scripts. Existing METIS
results and earlier experiments remain separate.

[upstream]: https://github.com/OpenQuantumComputing/min-cutrank/tree/86c1d8ff5d9b901983d9ad2254c8df0651e0219d
[archive]: ../archive/README.md#min-cutrank-comparison
[pdf]: results/min_cutrank/QuantumHamlets_min_cutrank_comparison.pdf
[docx]: results/min_cutrank/QuantumHamlets_min_cutrank_comparison.docx
[native]: results/min_cutrank/circuit_rank_native.png
[hybrid]: results/min_cutrank/hybrid_rank_30s.png
[evolution]: results/min_cutrank/evolution_rank.png
[remaining]: results/min_cutrank/remaining_circuits_native.png
[synthetic]: results/min_cutrank/synthetic_rank_native.png
[time]: results/min_cutrank/quality_vs_time.png
[sensitivity]: results/min_cutrank/schedule_comparison.json
[protocol]: results/min_cutrank/protocol.json
[summary]: results/min_cutrank/summary.tsv
[families]: results/min_cutrank/family_summary.tsv
[inputs]: results/min_cutrank/inputs.tsv
