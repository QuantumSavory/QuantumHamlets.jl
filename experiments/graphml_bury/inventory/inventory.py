#!/usr/bin/env python3
"""Inventory the supplied graphml* fixtures without running partition algorithms.

Only Python's standard library is used. XML is iterated, not loaded as a DOM.
Topology statistics describe the simple undirected projection; separate columns
record any directed edges, self loops, or duplicate edge declarations. Vertex
indices follow node declaration order. The 95th degree percentile uses nearest
rank (ceil(0.95*n)), with zero for an empty graph.
"""

import csv
import hashlib
import json
import math
import re
import statistics
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent


def tag(element):
    return element.tag.rsplit("}", 1)[-1]


def data_values(element):
    return {child.attrib["key"]: (child.text or "")
            for child in element if tag(child) == "data"}


def inspect(path):
    nodes = []
    declarations = []
    schemas = []
    node_values = {}
    edge_values = {}
    graph_defaults = []
    graph_count = 0
    explicit_directed = 0
    for _, element in ET.iterparse(path, events=("end",)):
        kind = tag(element)
        if kind == "key":
            schema = dict(element.attrib)
            schema["default"] = next((child.text for child in element
                                      if tag(child) == "default"), None)
            schemas.append(schema)
            element.clear()
        elif kind == "node":
            nodes.append(element.attrib["id"])
            for key, value in data_values(element).items():
                node_values.setdefault(key, Counter())[value] += 1
            element.clear()
        elif kind == "edge":
            declarations.append((element.attrib["source"], element.attrib["target"]))
            explicit_directed += element.attrib.get("directed", "false") in ("true", "1")
            for key, value in data_values(element).items():
                edge_values.setdefault(key, Counter())[value] += 1
            element.clear()
        elif kind == "graph":
            graph_defaults.append(element.attrib.get("edgedefault", ""))
            graph_count += 1
            element.clear()

    assert len(nodes) == len(set(nodes)), (path, "duplicate node IDs")
    assert graph_count == 1, (path, "expected one graph")
    mapping = {node: index for index, node in enumerate(nodes)}
    degrees = [0] * len(nodes)
    parents = list(range(len(nodes)))
    sizes = [1] * len(nodes)

    def find(v):
        while parents[v] != v:
            parents[v] = parents[parents[v]]
            v = parents[v]
        return v

    def union(u, v):
        u, v = find(u), find(v)
        if u != v:
            if sizes[u] < sizes[v]:
                u, v = v, u
            parents[v] = u
            sizes[u] += sizes[v]

    edges = set()
    loops = duplicates = 0
    for source, target in declarations:
        assert source in mapping and target in mapping, (path, source, target)
        u, v = mapping[source], mapping[target]
        if u == v:
            loops += 1
            continue
        edge = (min(u, v), max(u, v))
        if edge in edges:
            duplicates += 1
            continue
        edges.add(edge)
        degrees[u] += 1
        degrees[v] += 1
        union(u, v)

    n, m = len(nodes), len(edges)
    ordered_degrees = sorted(degrees)
    components = Counter(find(v) for v in range(n))
    numeric = all(re.fullmatch(r"-?\d+", node) for node in nodes)
    numeric_ids = list(map(int, nodes)) if numeric else []
    match = re.fullmatch(r"(.+?)_(\d+)", path.stem)
    gf_match = re.fullmatch(r"gf2\^(\d+)_mult", path.stem)
    family = match.group(1) if match else ("evolution" if path.parent.name == "graphml_evolution" else path.stem)
    parameter = int(match.group(2)) if match else (int(path.stem[1:]) if family == "evolution" else "")
    if gf_match:
        family, parameter = "gf2_mult", int(gf_match.group(1))
    topology = OUT / "topology" / (path.parent.name + "__" + path.stem + ".tsv")
    topology.parent.mkdir(exist_ok=True)
    with topology.open("w") as stream:
        stream.write(f"{n}\t{m}\n")
        for u, v in sorted(edges):
            stream.write(f"{u+1}\t{v+1}\n")
    node_map = OUT / "node_maps" / (path.parent.name + "__" + path.stem + ".tsv")
    node_map.parent.mkdir(exist_ok=True)
    with node_map.open("w", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(["vertex", "xml_id"])
        writer.writerows(enumerate(nodes, start=1))
    row = dict(
        directory=path.parent.name, graph=path.stem, family=family, parameter=parameter,
        path=str(path.relative_to(ROOT)), topology_path=str(topology.relative_to(ROOT)),
        node_map_path=str(node_map.relative_to(ROOT)), bytes=path.stat().st_size,
        n=n, m=m, declared_edges=len(declarations), density=2*m/(n*(n-1)) if n > 1 else 0,
        degree_min=min(degrees, default=0), degree_median=statistics.median(degrees) if n else 0,
        degree_mean=2*m/n if n else 0,
        degree_p95=ordered_degrees[math.ceil(0.95*n)-1] if n else 0,
        degree_max=max(degrees, default=0), isolated=degrees.count(0),
        components=len(components), largest_component=max(components.values(), default=0),
        self_loops=loops, duplicate_edges=duplicates,
        graph_default=graph_defaults[0], explicit_directed_edges=explicit_directed,
        numeric_ids=numeric,
        numeric_ids_ascending=(numeric_ids == sorted(numeric_ids)) if numeric else False,
        ids_are_1_to_n=numeric_ids == list(range(1, n+1)) if numeric else False,
        sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
    )
    meta = dict(path=row["path"], schemas=schemas,
                node_data={key: dict(counts) for key, counts in node_values.items()},
                edge_data={key: dict(counts) for key, counts in edge_values.items()},
                first_node_ids=nodes[:10], last_node_ids=nodes[-10:],
                minimum_numeric_id=min(numeric_ids, default=None),
                maximum_numeric_id=max(numeric_ids, default=None))
    return row, meta


def write_tsv(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main():
    paths = sorted((ROOT / "resources").glob("graphml*/*.graphml"))
    rows, metadata = [], []
    for path in paths:
        row, meta = inspect(path)
        rows.append(row)
        metadata.append(meta)
        print(f"{row['directory']}/{row['graph']}: n={row['n']} m={row['m']} "
              f"density={row['density']:.5f} mean_degree={row['degree_mean']:.2f} "
              f"max_degree={row['degree_max']} components={row['components']}", flush=True)
    write_tsv(OUT / "inventory.tsv", rows)
    (OUT / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    by_key = {(r['directory'], r['graph']): r for r in rows}
    pairs = []
    for initial in [r for r in rows if r['directory'] == 'graphml_initial']:
        reduced = by_key['graphml_reduced', initial['graph']]
        pairs.append(dict(graph=initial['graph'], n_initial=initial['n'], n_reduced=reduced['n'],
                          vertex_ratio=reduced['n']/initial['n'], m_initial=initial['m'],
                          m_reduced=reduced['m'], edge_ratio=reduced['m']/initial['m'],
                          density_initial=initial['density'], density_reduced=reduced['density'],
                          mean_degree_initial=initial['degree_mean'], mean_degree_reduced=reduced['degree_mean']))
    write_tsv(OUT / "initial_reduced.tsv", pairs)
    print(f"INVENTORY_COMPLETE files={len(rows)} total_bytes={sum(r['bytes'] for r in rows)}")


if __name__ == "__main__":
    main()
