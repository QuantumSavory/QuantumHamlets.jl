#!/usr/bin/env python3
"""Generate per-file property and backend guidance from audited results."""
import csv
from pathlib import Path

HERE = Path(__file__).resolve().parent
with (HERE/"inventory/inventory.tsv").open() as f:
    inventory = list(csv.DictReader(f, delimiter="\t"))
with (HERE/"comparison.tsv").open() as f:
    comparisons = {(r["directory"],r["graph"],r["kind"]):r
                   for r in csv.DictReader(f, delimiter="\t") if r["half"]=="True"}
out=["# Per-graph BURY backend guide", "",
     "All choices below concern the measured half-size target `floor(n/2)`. "
     "Seeded means explicit mapped vertex 1, the first declared XML node. "
     "Multiple backends mean their medians differ by at most 10% from the lowest; "
     "the first listed has the lowest median. This is a practical comparison band, "
     "not a statistical confidence interval. Seed or target changes can change the winner.", "",
     "All graphs are connected, simple, undirected, and have no isolated vertices. "
     "Density is `2m/(n(n-1))`. Complete numerical timings are in "
     "[comparison.tsv](comparison.tsv); additional degree properties are in "
     "[the inventory](inventory/inventory.tsv).", ""]
for directory in ("graphml_initial","graphml_reduced","graphml_evolution"):
    out += ["## "+directory, "",
            "| Graph | Vertices | Edges | Mean degree | Density | Ordinary | Explicit seeded |",
            "| --- | ---: | ---: | ---: | ---: | --- | --- |"]
    rs=sorted((r for r in inventory if r["directory"]==directory),
              key=lambda r:(r["family"],int(r["parameter"] or 0)))
    for r in rs:
        choices=[]
        for kind in ("ordinary","seeded"):
            c=comparisons.get((directory,r["graph"],kind))
            choices.append(c["within_10pct"].replace("/"," / ") if c else "not measured")
        out.append(f'| `{r["graph"]}` | {int(r["n"]):,} | {int(r["m"]):,} | '
                   f'{float(r["degree_mean"]):.2f} | {100*float(r["density"]):.3f}% | '
                   f'{choices[0]} | {choices[1]} |')
    out += [""]
out += ["Automatic search is separate: heap and scan are within 10% on 42 of 44 "
        "completed graphs. Automatic evolution `k5` stopped after a slow legacy "
        "preflight; no scan/heap automatic timing exists for it. See the "
        "[findings report](../../docs/graphml_backend_findings.md) for timing examples and limits.", ""]
(HERE/"per_graph_guide.md").write_text("\n".join(out))
print("Saved per_graph_guide.md for",len(inventory),"graphs")
