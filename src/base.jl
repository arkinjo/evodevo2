using Random
using DataStructures

@enum Face North East South West
@enum NAEnv Ancestral Novel
@enum TrainMode TestMode

const sixpii = 6.0/π
const sqrt3 = √3.0
function lcatan(x)
    sixpii*atan(x/sqrt3)
end

function relu(x)
    max(0, x)
end

mutable struct Setting
    basename::String
    seed::Int64
    with_cue::Bool
    max_pop::Int64
    num_cell_x::Int64
    num_cell_y::Int64
    num_env::Int64
    num_layers::Int64
    num_dev::Int64
    num_components::Vector{Int64}
    topology::Vector{OrderedDict{Int64,Float32}}
    omega::Vector{Float32}
    afuncs::Vector{Function}
    state_memory::Vector{Float32}
    env_noise::Float32
    mut_rate::Float32
    conv_dev::Float32
    denv::Float32
    selstrength::Float32
end

const density_default = 0.02
function default_setting(basename::String, num_layers::Int64=5, seed::Int64=13579)
    num_components = fill(200, num_layers)
    num_env = 200 ÷ 4
    topology = Vector()
    afuncs = Vector()
    for i = 1:num_layers
        push!(topology, OrderedDict())
        if i > 1 
            topology[i][i-1] = density_default
            if i < num_layers
                topology[i][i] = density_default
                push!(afuncs, lcatan)
            else
                push!(afuncs, tanh)
            end
        elseif i == 1
            push!(afuncs, identity)
        end

    end
    Setting(basename,
            seed,
            true, # with_cue
            200, # max_pop
            1,   # num_cell_x
            1,   # num_cell_y
            num_env,
            num_layers,   # num_layers
            200, # num_dev 
            num_components,      # num_components
            topology, 
            ones(Float32, num_layers), # omega (use "set_omegas")
            afuncs, 
            zeros(num_layers), # state_memory
            0.04, # env_noise
            0.001, # mut_rate
            1e-5, # conv_dev
            0.5,  # magnitude of environmental changes [0,1]
            20.0 # selstrength
            )
end

function set_omegas(s::Setting)
    for l = 2:s.num_layers
        omega = 0.0
        for (k,d) in s.topology[l]
            omega += d * s.num_components[k]
        end
        s.omega[l] = 1.0/sqrt(omega)
    end
    s
end
