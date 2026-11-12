# Vendored min-cutrank

The four algorithm files and MIT license are copied from
[OpenQuantumComputing/min-cutrank](https://github.com/OpenQuantumComputing/min-cutrank/tree/86c1d8ff5d9b901983d9ad2254c8df0651e0219d),
commit `86c1d8ff5d9b901983d9ad2254c8df0651e0219d` (the version used in the
comparison study). This folder is a Python package loaded by the Julia baseline
wrappers; no separate Python package installation is required.

Local changes to upstream files:

- Imports are package-relative so they do not replace unrelated Python modules.
- The two annealing functions accept an optional RNG. The Julia adapter supplies
  a local `random.Random(seed)` instead of modifying Python's global RNG.
- Scalar `numpy.exp` calls use standard-library `math.exp`. The core uses Python
  lists for its matrices; this removes its only NumPy dependency. Floating-point
  acceptance thresholds can differ at the last bit across math implementations.

`__init__.py` adds the small `anneal(adjacency, flags, temperatures, seed,
keep_best=True)` adapter. It copies inputs, checks nonempty capacities and
temperatures, invokes `cut_rank_annealing_row_formula`, and returns Boolean
membership. A subclass records the first strictly best completed state,
including the initializer. Set `keep_best=False` for the native final-state
policy. Recording candidates does not change proposals or RNG draws. Both part
sizes are preserved by swaps; the Julia wrapper chooses or supplies the initial
partition. There is no timeout or restart policy in this adapter.
An empty temperature sequence returns the initializer without running sweeps.

The graph/basis updates and swap-rank formulas remain unchanged. Copyright and
permission are retained in [LICENSE](LICENSE).
