# Completed research archives

Archived on 12 November 2026. These are complete preservation snapshots of the retired loose studies, with original relative paths. Active algorithms, the current Toffoli example, original graph collections, and personal research scripts remain in the working tree.

Each ZIP contains `ARCHIVE_MANIFEST.json` with file sizes, SHA-256 hashes, modes, original timestamps, and the captured package commit. All archives were extracted and every payload file compared with its source; a second local copy was also hash-checked. See [verification.json](verification.json) and [SHA256SUMS](SHA256SUMS).

The migration removed 3,325 verified loose payload files and 15 regenerable cache files, retaining 26 readable assets plus the required test reference. Post-migration package tests passed 459,283 assertions. A restored historical GraphML workflow and the maintained seven-method Toffoli example both loaded and ran their smoke checks; 92 local documentation links and 29 heading anchors were checked. Runtime source, original graph collections, and algorithm policies were preserved.

| Study | Archive | Payload files | ZIP size |
| --- | --- | ---: | ---: |
| [Bipartition range performance](#bipartition-range-performance) | [ZIP](bipartition_range_performance_20261112T212343Z.zip) | 8 | 0.10 MiB |
| [Bug audit](#bug-audit) | [ZIP](bug_audit_20261112T212343Z.zip) | 21 | 0.06 MiB |
| [BURY approximation](#bury-approximation) | [ZIP](bury_approximation_20261112T212343Z.zip) | 22 | 0.17 MiB |
| [BURY backends](#bury-backends) | [ZIP](bury_backends_20261112T212343Z.zip) | 9 | 0.03 MiB |
| [BURY performance](#bury-performance) | [ZIP](bury_performance_20261112T212343Z.zip) | 22 | 0.19 MiB |
| [Cost performance](#cost-performance) | [ZIP](cost_performance_20261112T212343Z.zip) | 12 | 0.02 MiB |
| [GraphML BURY](#graphml-bury) | [ZIP](graphml_bury_20261112T212343Z.zip) | 640 | 1.33 MiB |
| [Min-cutrank comparison](#min-cutrank-comparison) | [ZIP](min_cutrank_comparison_20261112T212343Z.zip) | 2,520 | 35.99 MiB |
| [GraphML comparisons](#graphml-comparisons) | [ZIP](graphml_comparisons_20261112T212343Z.zip) | 66 | 4.09 MiB |
| [Historical docs](#historical-docs) | [ZIP](historical_docs_20261112T212343Z.zip) | 5 | 0.03 MiB |
| [Restoration context](#restoration-context) | [ZIP](restoration_context_20261112T212343Z.zip) | 139 | 28.45 MiB |

The eleven ZIPs total 70.48 MiB. Counts include the shared restoration snapshot and nested upstream Git files, so they are not counts of unique research outputs. Caches (`__pycache__`, `.pyc`, `.CondaPkg`) and `.DS_Store` are explicitly excluded; excluded paths are recorded. Ignored PNG/PDF and untracked research files are included.

## Reading and restoring

For normal reading, use the maintained reports in `docs/` and their selected assets in `docs/results/`. Raw tasks, complete memberships, audits, scripts, and earlier figures are inside the relevant ZIP. Markdown links cannot open a member inside a ZIP.

For inspection, extract a study anywhere. For execution, extract its ZIP, its dependencies listed below, and `restoration_context` into the **same scratch package directory**. The internal paths begin with `experiments/`, `resources/`, `src/`, and so on. Avoid extracting over the active checkout. Work on a copy and keep the original archive unchanged.

Each ZIP has its own file named `ARCHIVE_MANIFEST.json`. Verify each archive separately before combining extracted trees, or retain these manifests under distinct study names. A combined extraction can overwrite that one metadata filename; the last manifest describes only its own ZIP, not the whole restored workspace. The original manifests always remain inside their respective ZIPs.

Check archive integrity before extraction (from this directory):

```sh
shasum -a 256 -c SHA256SUMS
```

If a historical script requires package Git metadata, first extract the context to a temporary folder, clone its `archival_context/repository.bundle` into a separate scratch checkout, then overlay the required ZIPs there. For a particular historical experiment, use its recorded commit/snapshot and environment rather than assuming the archival-time working tree is the right version. Original run logs, hashes, and snapshots are preserved without retroactively rewriting provenance.

The parent environment snapshot references other local checkouts (including FBQC and QuantumClifford/QECCore); it is not a portable standalone environment. Restore those dependencies or build a separate research environment from the recorded requirements. The old benchmark manifest also predates the new SA dependencies. Do not resolve or edit the only preserved copy of a historical manifest. Report regeneration may require the original external document-rendering tools; saved PDF/DOCX files can be read directly. Archiving preserves available evidence and code, not a guarantee of indefinite third-party runtime availability.

For file-level verification after extraction, read `ARCHIVE_MANIFEST.json` and compare the SHA-256 and byte length of each listed path. Original modes and nanosecond timestamps are recorded there because ZIP tools do not uniformly restore them. The preserved one-time packaging script documents the checks performed.

## Local recovery copy

A second verified local copy is stored at `/Users/micciche/Research/QuantumInfo23/JuliaProjects/QuantumHamlets-archive-backup/20261112T212343Z`. It is on the same computer, not an off-device backup. The repository ZIPs are the portable copies; keep them with a normal repository/data backup. These archives reduce loose-file clutter without rewriting Git history. Previously committed files remain in Git history.

## Bipartition range performance

The target-range benchmark, frozen old wrapper, results, and full package validation log. Restore `bury_performance` for its benchmark environment.

## Bug audit

Original bug/cost audits, snapshots, fixtures, and logs. Restore `bury_performance` for its environment. Some historical audit scripts deliberately assert old failures or overwrite their output; run copies in a scratch checkout.

## BURY approximation

Constructors, proofs-related checks, certified/search prototypes, logs, and original figures. Read the maintained [approximation discussion](../docs/bury_approximation.md).

## BURY backends

Integrated-backend checks and timing evidence. Restore `bury_performance` and `graphml_bury` as well: its driver reads their scripts and inventory.

## BURY performance

Prototype implementations, benchmark environment, raw measurements, profiles, and validations. Shared context includes the two required QAOA fixtures. Read the [findings](../docs/bury_performance_findings.md).

## Cost performance

Frozen evaluator reference, before/after results, and scripts. Restore `bury_performance` for the benchmark environment. The active package test now uses an identical copy at `test/reference/costs_before_optimization.jl`.

## GraphML BURY

All inventory/topology exports, prototype jobs, samples, logs, and backend comparisons. Restore `bury_performance` for its implementations/environment. Read the [findings](../docs/graphml_backend_findings.md) and [per-graph guide](../docs/results/graphml_bury/per_graph_guide.md).

## Min-cutrank comparison

The complete study, including all task JSON files, raw scores and memberships, schedule sensitivity, tiny fixtures, original and final figures, report drafts/QA, both upstream tarballs, the earlier curated ZIP, and the nested upstream `.git` history. Restore `graphml_bury` and `graphml_comparisons` too: input export reads their inventories and saved results. Read the [summary](../docs/min_cutrank_comparison.md) or [full PDF](../docs/results/min_cutrank/QuantumHamlets_min_cutrank_comparison.pdf). The original run is `run_20261112T161827Z`; the newer callable Julia baselines and Toffoli example remain active.

## GraphML comparisons

Both historical BURY/METIS GraphML workflows and all their result directories under `resources/figs/graphml_comparison/` and `graphml_families/`, including stopped and superseded runs. Both script families are together because the newer analysis includes the older module. The original guide is `resources/graphml_family_plots_README.md` inside the ZIP. Their native METIS policy used 3% imbalance; archival does not change the findings or balance interpretation.

## Historical docs

Five original documents: `bipartition_range_handoff.md`, `bipartition_range_review.md`, `bug_audit.md`, `library_experiments_handoff.md`, and `revision_proposal.md`. Their original paths are retained under `docs/`. Current limitations remain in the [usage guide](../docs/usage.md); remaining cleanup proposals are in the [streamlining plan](../docs/streamlining_plan.md).

## Restoration context

The archival-time `src/`, `test/`, package environment, Markdown documentation, resource scripts, 45 circuit GraphML inputs, two benchmark QAOA inputs, and top-level example figures. `archival_context/` contains a verified Git bundle with all available package refs, the pre-archive dirty-tree status/diff, the one-time packaging script, and copies of the current parent Project/Manifest. This is dated context, not a claim that today's source or environment was used for every old run.
