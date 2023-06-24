include("base.jl")
include("env.jl")
include("genome.jl")
include("cell.jl")
include("indiv.jl")
include("population.jl")

s = default_setting
if ARGS[1] == "nocue"
    s.with_cue = false
end
s = set_omegas(s)
println("#with_cue= $(s.with_cue)")

train_epochs(10, 200, s)
