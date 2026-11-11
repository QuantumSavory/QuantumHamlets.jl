using Test
using Graphs
using QuantumHamlets
using Random

const Q = QuantumHamlets

include("legacy.jl")
include("known_bugs.jl")
include("api.jl")
include("bury_backends.jl")
include("bury_validation.jl")
include("baseline_legacy_regressions.jl")
include("cost_regressions.jl")
include("research_scores.jl")
include("cost_semantics_regressions.jl")
include("rank_regressions.jl")
