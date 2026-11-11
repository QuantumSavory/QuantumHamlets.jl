"""Check the matching-cut/component-packing characterization on tiny graphs.

Standalone mathematical prototype, not a production BURY backend.  An exact
small-separator oracle costs sum_{j=0}^q binomial(n,j) component/DP passes.
Run: python3 experiments/bury_approximation/theory/separator_oracle.py
"""
from itertools import combinations
import random


def matching_cost(adj, side):
    """Independent augmenting-path bipartite matching."""
    mate = {}
    def augment(u, seen):
        for v in adj[u] - side:
            if v in seen:
                continue
            seen.add(v)
            if v not in mate or augment(mate[v], seen):
                mate[v] = u
                return True
        return False
    return sum(augment(u, set()) for u in sorted(side))


def pack_components(adj, removed, target):
    """Return a target-sized side whose cut is covered by removed, if any."""
    unseen = set(range(len(adj))) - removed
    components = []
    while unseen:
        root = min(unseen)
        unseen.remove(root)
        component, stack = {root}, [root]
        while stack:
            u = stack.pop()
            fresh = adj[u] & unseen
            unseen -= fresh
            component |= fresh
            stack.extend(fresh)
        components.append(component)
    attainable = {0: set()}
    for component in components:
        for size, union in list(attainable.items()):
            new_size = size + len(component)
            if new_size <= target and new_size not in attainable:
                attainable[new_size] = union | component
    for size in range(max(0, target - len(removed)), target + 1):
        if size in attainable:
            return attainable[size] | set(sorted(removed)[:target-size])
    return None


def separator_oracle(adj, target, q):
    for size in range(q + 1):
        for vertices in combinations(range(len(adj)), size):
            removed = set(vertices)
            side = pack_components(adj, removed, target)
            if side is not None:
                return side, removed
    return None


def check_graph(adj):
    n = len(adj)
    assertions = 0
    for target in range(n + 1):
        optimum = min(matching_cost(adj, set(vs))
                      for vs in combinations(range(n), target))
        for q in range(min(target, n-target) + 1):
            result = separator_oracle(adj, target, q)
            assert (result is not None) == (optimum <= q)
            assertions += 1
            if result is not None:
                side, removed = result
                assert len(side) == target
                assert len(removed) <= q
                assert matching_cost(adj, side) <= len(removed)
                assert all(u in removed or v in removed
                           for u in side for v in adj[u] - side)
                assertions += 4
    return assertions


def main():
    graphs = assertions = 0
    for n in range(6):
        pairs = list(combinations(range(n), 2))
        for bits in range(1 << len(pairs)):
            adj = [set() for _ in range(n)]
            for i, (u, v) in enumerate(pairs):
                if bits >> i & 1:
                    adj[u].add(v)
                    adj[v].add(u)
            assertions += check_graph(adj)
            graphs += 1
    rng = random.Random(20261111)
    for n in range(6, 11):
        for _ in range(20):
            p = rng.random()
            adj = [set() for _ in range(n)]
            for u, v in combinations(range(n), 2):
                if rng.random() < p:
                    adj[u].add(v)
                    adj[v].add(u)
            assertions += check_graph(adj)
            graphs += 1
    print(f"PASS: {graphs} graphs; {assertions} assertions")


if __name__ == "__main__":
    main()
