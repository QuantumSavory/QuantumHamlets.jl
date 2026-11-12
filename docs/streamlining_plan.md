# Proposal to streamline QuantumHamlets.jl

**12 November 2026 — archival phase complete; broader cleanup remains proposed**

The aim is a library that is easy to read, call, and extend, with a small set of
research examples. Start by archiving completed research work and consolidating
documentation. Then coalesce small source and test files while preserving
current behavior. Decisions to retire algorithms or APIs should be separate.

This report records the recommendations discussed with the user. The user
subsequently authorized the archival phase. Completed snapshots and restoration
instructions are listed in the [archive index](../archive/README.md).
The completed phase removed verified loose study copies, retained readable
deliverables, repaired links, and moved the required cost reference into tests.
Source coalescing, dependency removal, and API retirement remain proposals.

## 1. Where the clutter is

The working-tree inventory, excluding Git internals and Python/environment
caches, is approximately:

| Area | Files | Size | Interpretation |
| --- | ---: | ---: | --- |
| `src/` | 21 | 0.1 MiB | Small runtime implementation, including vendored Python source. |
| `test/` | 16 | 0.1 MiB | Useful coverage, with several small historical regression files. |
| `docs/` | 28 | 12.9 MiB | Guides, research reports, old planning documents, and figures. |
| `experiments/` | 3,228 | 107.9 MiB | Main source of file-count growth. |
| `resources/` | 219 | 160.3 MiB | Graph inputs, examples, historical workflows, and generated figures. |

These counts precede this report and describe files on disk, not only tracked
files. They are a prioritization aid rather than a future target.

The completed min-cutrank comparison accounts for **2,494 files / 90.8 MiB**,
including 2,010 individual annealing task JSON files. The original QAOA graph
collections occupy about 144.8 MiB but contain only 83 files. Archiving completed
experiments therefore addresses the file-count problem more directly than
moving the original datasets.

## 2. Recommended maintained surface

Keep these easy to find:

- The package source and its tests.
- `README.md` and [the usage guide](usage.md).
- [The BURY approximation discussion](bury_approximation.md) and
  [the min-cutrank comparison report](min_cutrank_comparison.md).
- The short [Toffoli comparison function](../resources/toffoli_baselines.jl)
  and useful grid/QAOA examples.
- Original graph datasets and a small selection of final figures and tables
  referenced by maintained reports.

Keep the five convenient baseline calls: `min_cutrank`, `bury_sa`,
`bury_seeded_sa`, `metis_partition`, and `metis_sa`. A larger dispatcher would
save few lines and make research scripts less readable. The intended pattern
remains a direct function call returning `GraphPartition`.

Do not replace the existing experiment machinery with another generic experiment
framework. A few readable Julia examples should be sufficient for new work.

## 3. Source consolidation

These changes reduce top-level Julia source files from ten to seven without
changing function names, signatures, return values, or partition policies:

| Current file | Proposed destination | Reason |
| --- | --- | --- |
| `src/display.jl` | `src/quantum_land.jl` | Its 14 lines display those legacy objects. |
| `src/graph_utils.jl` | BURY helpers into `src/bury.jl`; legacy edge removal beside legacy objects | Three helpers do not need a miscellaneous file. |
| `src/api.jl` | Partitioning functions into `src/bury.jl`, cost adapters into `src/costs.jl`, plotting adapter into `src/visualization.jl` | Put each public function beside its implementation. |

Retain `src/partition.jl` for `GraphPartition` and its accessors. Keep heap and
scan implementations separate; they are substantial algorithms rather than
scaffolding. Preserve the structure of `src/min_cutrank/`, its upstream version
record, local-change notes, and MIT license. Flattening its four upstream
modules would make future comparison with upstream harder for little benefit.

## 4. Dependencies and possible feature retirement

| Candidate | Recommendation | Prerequisite or tradeoff |
| --- | --- | --- |
| Nemo | Move to test-only dependencies. | Production rank uses packed GF(2) elimination. Reference tests must import Nemo directly instead of through QuantumHamlets. |
| Colors | Remove its redundant direct import/dependency. | No direct production use was found; plotting dependencies retain their own color support. |
| Saran–Vazirani | Consider retiring if no longer used. | It alone uses GraphsFlows, MultiwayNumberPartitioning, and HiGHS. Removal changes a supported baseline. |
| Runtime `backend=:legacy` | Consider retiring in a later API decision. | Documented calls and historical scripts use it. Preserve the frozen original under `test/reference/` regardless. |
| `quantumLand`, `villager`, `hamlet` | Preserve in the first pass. | Existing API, seed search, plotting, and references still use them. A migration is more than deleting a file. |
| `antipode_seeds` | Candidate for removal after preserving its source. | No current callers were found outside the frozen historical implementation; it is not in the public list. |
| IGraphs | Keep. | Current matching kernels depend on it. |

Keep exhaustive `best_seed` unless explicitly retired; it remains a public,
tested feature even though recent plots use sqrt sampling.

Making plotting or PythonCall optional could reduce installation and startup
weight, but would add extension files and usage rules. Defer that change unless
startup or dependency weight becomes a separate priority. Fewer files and fewer
mandatory dependencies are related but distinct goals.

## 5. Documentation, examples, and tests

Archive historical handoff prompts, the old revision proposal, implementation
reviews, and resolved bug audits. First carry any still-current limitations or
required explanations into the usage guide. The two backend performance reports
could become one concise report, with their originals preserved in the archive.

The older `graphml_plots*` / `run_graphml_plots.py` workflow and the later
`graphml_family_*` / `run_graphml_family_plots.py` workflow should become frozen
historical studies unless their broader multipartition analyses remain active.
The newer workflow includes the older one for shared helpers: archive them
together, or separate those helpers before retiring either workflow.

Keep `resources/toffoli_baselines.jl` as the current short example. Bring useful
grid/QAOA examples to the same pattern: a callable function, no experiment on
include, explicit arguments, and returned data. Review `sandbox.jl` and other
personal research scripts with the user before treating them as disposable.

Consolidate tests by responsibility while preserving assertions:

- Fold the small `known_bugs.jl` cases into the relevant API, baseline, or cost
  tests; their historical label no longer helps navigation.
- Combine overlapping cost regression files.
- Keep BURY backend equivalence, target-range behavior, and baseline tests
  identifiable rather than building one enormous test file.
- Preserve `test/reference/bury_original.jl` unchanged.

One dependency must be fixed before archiving experiments:
`test/cost_semantics_regressions.jl` includes
`experiments/cost_performance/reference_costs.jl`. Copy that frozen reference
into `test/reference/` without changing its contents, update the include, and
make its test-only imports explicit.

## 6. What archiving would entail

Archiving turns a completed study into a preserved snapshot that can be opened
when needed. It removes thousands of loose files from the everyday working
tree while keeping the evidence and reproduction material together.

### Archive contents

Use one archive per coherent completed study. Preserve the original relative
paths inside it, including files spread across `experiments/`, `resources/`,
and `docs/` when they belong to that study. Include:

- Raw results, individual tasks, memberships, inputs, seeds, schedules, and
  balance/timeout policies.
- Analysis and plotting code, reference implementations, upstream source and
  licenses, available source snapshots, and recorded version/environment data.
- Final and superseded figures, reports, audits, logs, preliminary/stopped runs,
  and saved inventories. These retain the history of the conclusions.
- Relevant untracked files and ignored PNG/PDF files. A tracked-files-only
  `git archive` would miss part of the current evidence.

Exclude regenerable caches, `.CondaPkg` environments, `__pycache__`, and operating
system metadata from the research archive; record the exclusion rule. Treat
existing nested repositories explicitly: retain required source/history rather
than silently losing it through a generic file filter.

Preserve the source/environment records captured for the original run. A
snapshot of today's package or parent environment is useful context but does
not replace the version actually used for the historical experiment. Record
missing provenance honestly. Preserving files does not guarantee that all
third-party environments will remain runnable indefinitely.

### Storage and everyday access

The simplest first arrangement is a small `archive/` directory containing one
ZIP per study and one index. ZIP is convenient for ordinary file browsing;
the archive format itself is not an important architectural decision.

For example, the index could identify a study by date, describe its scope,
link its readable report, name its archive and checksum, and give its original
paths and restoration instructions. Keep the current Markdown reports and
their selected figures/tables outside the archive so normal reading does not
require extraction. Their small duplication inside the complete snapshot is
intentional. Update report links to the retained files or archive index;
ordinary Markdown links cannot browse files inside a ZIP.

An archive within the repository reduces visible file count but not necessarily
Git storage. Previously committed files remain in Git history, and committing
a ZIP adds another representation. If smaller clones become a goal, store the
archive in a separately backed-up location with a stable retrieval link. A
sibling folder is convenient locally but does not travel with a clone. Do not
rewrite Git history as part of this cleanup.

### Verification before removing loose files

1. **Inventory the scope.** Identify all related files and external dependencies,
   including references from active tests, scripts, and documentation. Record
   pending changes so the current work is not mistaken for the last commit.
2. **Build a new snapshot.** Write to a fresh archive outside the directories
   being collected. Retain paths, file contents, and relevant file metadata.
   Do not overwrite the existing deliverable bundle.
3. **Record integrity information.** Store an internal manifest with each file's
   path, size, and SHA-256 digest, plus an external digest of the archive itself.
4. **Test restoration.** Extract into a temporary directory and compare every
   file with the manifest and source inventory. Check expected counts and
   inspect representative tables, plots, reports, and source files.
5. **Check usability.** Ensure the restoration instructions locate inputs,
   code, and environments. Run a small loading/analysis check where practical;
   a full expensive experiment rerun is unnecessary for validating a copy.
6. **Prepare the working-tree change.** Retain chosen readable deliverables,
   relocate required test references, update links, and list the exact loose
   files proposed for removal. Verify a separately recoverable archive copy.
7. **Review and apply the cleanup.** Remove loose copies only after the snapshot
   and removal scope are accepted. Run package tests and maintained examples,
   then check links again.

The existing min-cutrank artifact ZIP is a curated sharing bundle, **not a full
archive of its working directory**. Its packaging code intentionally excludes
superseded `figures/`, internal report-render QA and drafts, and the full
upstream-source tarball. It is useful, but cannot alone justify deleting all
loose study files.

### Returning to an archived study

Extract into a separate research checkout or scratch folder, following the
recorded layout. Work on a copy and give subsequent results a new study/run
identity. The original archive remains immutable. Keep the files loose again
if that study becomes active; archiving need not make future work cumbersome.

## 7. Proposed order and acceptance criteria

1. **Pilot the min-cutrank archive.** It offers the largest file-count reduction
   and has a clear current report and a replacement callable baseline API.
2. **Archive the other completed studies and historical workflows.** Preserve
   cross-directory dependencies and keep original graph collections available.
3. **Simplify documentation and source layout.** Apply the three source merges
   and reorganize tests without changing behavior.
4. **Remove unused runtime dependencies.** Move Nemo to tests and remove the
   redundant Colors dependency, then run the package suite.
5. **Decide on feature retirement separately.** Resolve Saran–Vazirani, runtime
   legacy BURY, and any longer-term legacy-object migration with the user.

The first pass succeeds when the everyday tree is easier to navigate, current
calls still work, maintained examples run, reports have working links, tests
pass, and every archived result can be restored with verified contents.
It should not change seeded search, target-range ties, cut definitions, METIS
balance policy, or annealing output policy as a side effect of reorganization.
