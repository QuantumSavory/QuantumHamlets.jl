# Streamlining plan and implementation record

**12 November 2026 — archiving complete; source consolidation implemented**

The aim is a library that is easy to read, call, and extend, with a small set of
research examples. Completed research has been archived and three source files
have been coalesced while preserving the existing API. Test regrouping,
dependency changes, and decisions to retire algorithms or APIs remain separate.

This report records the recommendations discussed with the user. The user
subsequently authorized archiving and then source consolidation. Snapshots and restoration
instructions are listed in the [archive index](../archive/README.md).

## Completed changes

**Archiving:** eight completed experiment collections, the coupled historical
GraphML workflows, five historical documents, and shared restoration context
were preserved in **11 ZIPs totaling 70.48 MiB**. All **3,464 archived payload
copies** were extracted and hash-checked, and a second local archive copy was
verified. These counts include shared context and nested Git files, rather than
3,464 independent research outputs.

The migration removed **3,325 verified loose payload files** and 15 regenerable
cache files. It retained 26 readable reports, plots, and summary files under
`docs/results/`, updated links, and copied the required cost reference
byte-for-byte into `test/reference/costs_before_optimization.jl`. The test now
uses that copy and imports its Nemo reference operations directly. Scoped
`.gitignore` exceptions allow the retained PNG/PDF deliverables to be tracked.
Runtime source, original graph collections, and algorithm policies were
preserved during this phase.

Post-archive validation passed **459,283 package assertions**. A restored
historical GraphML workflow and the maintained seven-method Toffoli example
passed smoke checks; 92 local documentation links and 29 anchors were checked.
Details and backup location are in the
[archival verification record](../archive/verification.json).

**Source consolidation:** the three moves in section 3 reduce top-level Julia
files from **ten to seven**, and total Julia source files from **thirteen to
ten**. Public functions now live beside their implementations. The module
include order loads legacy objects, `GraphPartition`, costs, BURY, baselines,
and plotting in that order. Function bodies and docstrings are moved with their
functions; exports, signatures, policies, dependencies, backend implementations,
vendored Python code, and frozen test references are preserved. The usage guide
now documents the actual file layout.

Post-consolidation validation passed **459,283 package assertions**, using:

```sh
julia --startup-file=no --project=. -e 'using Pkg; Pkg.test(; allow_reresolve=false)'
```

An independent parsed-source comparison confirmed that all **114 top-level
code and docstring expressions** were preserved exactly, apart from source
location metadata. The module's only change is its include list/order.
The edited documentation's 18 local links and heading anchors were checked.
These checks validate this consolidation separately from the archival work.

## 1. Where the clutter is

The **pre-cleanup** working-tree inventory, excluding Git internals and
Python/environment caches, was approximately:

| Area | Files | Size | Interpretation |
| --- | ---: | ---: | --- |
| `src/` | 21 | 0.1 MiB | Small runtime implementation, including vendored Python source. |
| `test/` | 16 | 0.1 MiB | Useful coverage, with several small historical regression files. |
| `docs/` | 28 | 12.9 MiB | Guides, research reports, old planning documents, and figures. |
| `experiments/` | 3,228 | 107.9 MiB | Main source of file-count growth. |
| `resources/` | 219 | 160.3 MiB | Graph inputs, examples, historical workflows, and generated figures. |

These historical counts precede the changes recorded above and describe files
on disk, not only tracked files. They explain the original priorities and are
not a current inventory or a future target.

The completed min-cutrank comparison accounted for **2,494 files / 90.8 MiB**,
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

The following moves are implemented. The three former files and their module
includes have been removed; their functions remain available under the same
names, signatures, return values, and partition policies:

| Former file | Current destination | What moved |
| --- | --- | --- |
| `src/display.jl` | `src/quantum_land.jl` | The three `Base.show` methods for legacy objects. |
| `src/graph_utils.jl` | `src/bury.jl` and `src/quantum_land.jl` | `_check_partition_count` and `inverse_map` accompany BURY; `remove_local_edges` accompanies legacy objects. |
| `src/api.jl` | `src/bury.jl`, `src/costs.jl`, `src/partition.jl`, and `src/visualization.jl` | Partition functions and target helpers accompany BURY; scalar cost adapters and `matching_details` accompany cost kernels; `_labels_for` accompanies `GraphPartition`; `_land_for` and `plot_partition` accompany plotting. |

`src/partition.jl` retains `GraphPartition`, its conversions, and accessors.
Heap, scan, and legacy implementations remain separate under `src/bury/`.
The structure of `src/min_cutrank/`, upstream version record, local-change notes,
and MIT license are unchanged. Existing package calls need no migration.
Historical archive snapshots keep their original source layout for restoration.

## 4. Dependencies and possible feature retirement

These are remaining proposals; the source consolidation does not remove or
upgrade dependencies or retire any algorithm/API.

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

Historical handoff prompts, the old revision proposal, implementation reviews,
and resolved bug audits are now in the historical documentation archive. The
usage guide retains current limitations. Combining the two backend performance
reports into one concise report remains a possible follow-up; both are still
directly readable and their originals are preserved in the archives.

The older `graphml_plots*` / `run_graphml_plots.py` workflow and the later
`graphml_family_*` / `run_graphml_family_plots.py` workflow are frozen together
in the GraphML comparisons archive. The newer workflow includes the older one
for shared helpers; their original relative layout and result directories are
preserved together.

Keep `resources/toffoli_baselines.jl` as the current short example. Bring useful
grid/QAOA examples to the same pattern: a callable function, no experiment on
include, explicit arguments, and returned data. Review `sandbox.jl` and other
personal research scripts with the user before treating them as disposable.

Test regrouping has not been performed. Possible follow-up work, preserving
assertions, is to:

- Fold the small `known_bugs.jl` cases into the relevant API, baseline, or cost
  tests; their historical label no longer helps navigation.
- Combine overlapping cost regression files.
- Keep BURY backend equivalence, target-range behavior, and baseline tests
  identifiable rather than building one enormous test file.
- Preserve `test/reference/bury_original.jl` unchanged.

The archive dependency has been fixed: `test/cost_semantics_regressions.jl` now
includes `test/reference/costs_before_optimization.jl`. That file is an unchanged
copy of the former `experiments/cost_performance/reference_costs.jl`. Its
reference imports are explicit, and the full package suite passed after the
relocation. The frozen BURY reference remains unchanged.

## 6. Archive arrangement and restoration policy

Archiving turns a completed study into a preserved snapshot that can be opened
when needed. It removes thousands of loose files from the everyday working
tree while keeping the evidence and reproduction material together. This phase
is complete. The process below records the preservation policy used and the
policy to follow for any future archived studies.

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

The chosen arrangement is `archive/`, containing the study ZIPs, shared-context
and historical-document ZIPs, an index, checksums, and a verification record.
ZIP is convenient for ordinary file browsing; the format itself is not an
important architectural decision.

The index identifies the studies and dates, describes their scope, links readable
reports, names the archives and checksums, and gives restoration instructions.
Current Markdown reports and selected figures/tables remain outside the ZIPs so
normal reading does not require extraction. Their small duplication inside the
complete snapshots is intentional. Report links now target the retained files
or archive index; ordinary Markdown links cannot browse files inside a ZIP.

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

The earlier min-cutrank artifact ZIP was a curated sharing bundle, **not a full
archive of its working directory**. Its packaging code intentionally excludes
superseded `figures/`, internal report-render QA and drafts, and the full
upstream-source tarball. The new complete min-cutrank archive preserves that
bundle along with the omitted material, including the nested upstream Git
repository. The curated bundle was not used as the sole preservation copy.

### Returning to an archived study

Extract into a separate research checkout or scratch folder, following the
recorded layout. Work on a copy and give subsequent results a new study/run
identity. The original archive remains immutable. Keep the files loose again
if that study becomes active; archiving need not make future work cumbersome.

## 7. Current status and remaining work

| Work | Status |
| --- | --- |
| Archive completed studies, workflows, and historical documentation | Complete; verified snapshots, retained deliverables, repaired links, and restoration instructions are available. |
| Relocate the test reference needed after archiving | Complete; copied unchanged and validated by the full package suite. |
| Coalesce the three source files and update the usage guide | Complete; 459,283 assertions pass and parsed-source comparison confirms preserved definitions/docstrings. |
| Regroup test files and combine performance reports | Proposed; not part of this source-only consolidation. |
| Move Nemo to tests and remove the redundant Colors dependency | Proposed; package dependencies remain unchanged. |
| Retire Saran–Vazirani, runtime legacy BURY, or legacy objects | Separate decisions requiring discussion with the user. |
| Optional plotting or PythonCall loading | Deferred unless startup/dependency weight becomes a priority. |

The first pass succeeds when the everyday tree is easier to navigate, current
calls still work, maintained examples run, reports have working links, tests
pass, and every archived result can be restored with verified contents.
It should not change seeded search, target-range ties, cut definitions, METIS
balance policy, or annealing output policy as a side effect of reorganization.
