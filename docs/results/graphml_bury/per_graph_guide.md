# Per-graph BURY backend guide

All choices below concern the measured half-size target `floor(n/2)`. Seeded means explicit mapped vertex 1, the first declared XML node. Multiple backends mean their medians differ by at most 10% from the lowest; the first listed has the lowest median. This is a practical comparison band, not a statistical confidence interval. Seed or target changes can change the winner.

All graphs are connected, simple, undirected, and have no isolated vertices. Density is `2m/(n(n-1))`. Complete numerical timings are in [comparison.tsv](comparison.tsv); additional degree properties are in [the inventory](inventory.tsv).

## graphml_initial

| Graph | Vertices | Edges | Mean degree | Density | Ordinary | Explicit seeded |
| --- | ---: | ---: | ---: | ---: | --- | --- |
| `barenco_tof_3` | 83 | 102 | 2.46 | 2.997% | heap | scan |
| `barenco_tof_4` | 151 | 192 | 2.54 | 1.695% | heap | heap / scan |
| `barenco_tof_5` | 221 | 284 | 2.57 | 1.168% | heap | heap |
| `barenco_tof_10` | 571 | 744 | 2.61 | 0.457% | heap | heap |
| `csla_mux_3` | 239 | 296 | 2.48 | 1.041% | heap | heap |
| `csum_mux_9` | 622 | 760 | 2.44 | 0.394% | heap | heap |
| `gf2^4_mult` | 306 | 393 | 2.57 | 0.842% | heap | heap |
| `gf2^5_mult` | 465 | 604 | 2.60 | 0.560% | heap | heap |
| `gf2^6_mult` | 652 | 855 | 2.62 | 0.403% | heap | heap |
| `hwb6` | 325 | 434 | 2.67 | 0.824% | heap | heap |
| `mod_mult_55` | 167 | 206 | 2.47 | 1.486% | heap | heap |
| `mod_red_21` | 355 | 449 | 2.53 | 0.715% | heap | heap |
| `qcla_com_7` | 605 | 755 | 2.50 | 0.413% | heap | heap |
| `qft_4` | 143 | 184 | 2.57 | 1.812% | heap | heap / scan |
| `rc_adder_6` | 262 | 331 | 2.53 | 0.968% | heap | heap |
| `tof_3` | 67 | 80 | 2.39 | 3.618% | heap | scan |
| `tof_4` | 107 | 130 | 2.43 | 2.292% | heap | scan |
| `tof_5` | 147 | 180 | 2.45 | 1.677% | heap | heap |
| `tof_10` | 347 | 430 | 2.48 | 0.716% | heap | heap |
| `vbe_adder_3` | 200 | 248 | 2.48 | 1.246% | heap | heap |

## graphml_reduced

| Graph | Vertices | Edges | Mean degree | Density | Ordinary | Explicit seeded |
| --- | ---: | ---: | ---: | ---: | --- | --- |
| `barenco_tof_3` | 25 | 50 | 4.00 | 16.667% | scan | scan |
| `barenco_tof_4` | 39 | 91 | 4.67 | 12.281% | scan | scan |
| `barenco_tof_5` | 53 | 132 | 4.98 | 9.579% | heap / scan | scan |
| `barenco_tof_10` | 123 | 341 | 5.54 | 4.545% | heap | scan / heap |
| `csla_mux_3` | 111 | 237 | 4.27 | 3.882% | heap | scan |
| `csum_mux_9` | 166 | 400 | 4.82 | 2.921% | heap | heap / scan |
| `gf2^4_mult` | 132 | 415 | 6.29 | 4.800% | heap | scan / heap |
| `gf2^5_mult` | 225 | 882 | 7.84 | 3.500% | heap | scan / heap |
| `gf2^6_mult` | 294 | 1,186 | 8.07 | 2.754% | heap | heap |
| `hwb6` | 137 | 480 | 7.01 | 5.152% | heap | scan |
| `mod_mult_55` | 73 | 195 | 5.34 | 7.420% | heap | scan |
| `mod_red_21` | 133 | 353 | 5.31 | 4.021% | heap | scan / heap |
| `qcla_com_7` | 169 | 402 | 4.76 | 2.832% | heap | heap |
| `qft_4` | 87 | 158 | 3.63 | 4.223% | heap | scan |
| `rc_adder_6` | 106 | 226 | 4.26 | 4.061% | heap | scan |
| `tof_3` | 27 | 48 | 3.56 | 13.675% | scan | scan |
| `tof_4` | 41 | 72 | 3.51 | 8.780% | scan / heap | scan |
| `tof_5` | 55 | 96 | 3.49 | 6.465% | heap | scan |
| `tof_10` | 125 | 241 | 3.86 | 3.110% | heap | heap |
| `vbe_adder_3` | 44 | 98 | 4.45 | 10.359% | scan / heap | scan |

## graphml_evolution

| Graph | Vertices | Edges | Mean degree | Density | Ordinary | Explicit seeded |
| --- | ---: | ---: | ---: | ---: | --- | --- |
| `k1` | 25 | 42 | 3.36 | 14.000% | scan | scan |
| `k2` | 159 | 498 | 6.26 | 3.965% | heap | scan |
| `k3` | 562 | 2,166 | 7.71 | 1.374% | heap | heap |
| `k4` | 1,159 | 6,216 | 10.73 | 0.926% | heap | heap |
| `k5` | 1,772 | 12,102 | 13.66 | 0.771% | heap | heap |

Automatic search is separate: heap and scan are within 10% on 42 of 44 completed graphs. Automatic evolution `k5` stopped after a slow legacy preflight; no scan/heap automatic timing exists for it. See the [findings report](../../graphml_backend_findings.md) for timing examples and limits.
