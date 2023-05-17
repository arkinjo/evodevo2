module Hoge
include("base.jl")
include("genome.jl")

ds = default_setting
genome1 = make_genome(ds.num_layers, ds.num_components, ds.feedback_loops)
genome2 = make_genome(ds.num_layers, ds.num_components, ds.feedback_loops)
jldopen("genome.jld2","w") do file
    file["genome1"] = genome1
end
jldopen("genome.jld2","a+") do file
    file["genome2"] = genome2
end

end
