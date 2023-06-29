include("base.jl")
include("env.jl")
include("genome.jl")
include("cell.jl")
include("indiv.jl")
include("population.jl")

s = default_setting()
if ARGS[1] == "nocue"
    s.with_cue = false
end
s = set_omegas(s)
s.topology[2][2] = 0.02
s.topology[3][3] = 0.02
s.topology[4][4] = 0.02
s.state_memory[2] = 0.0
s.max_pop = 200
println("#with_cue= $(s.with_cue)")
println("#topology= $(s.topology)")

train_epochs(10, 200, s)
