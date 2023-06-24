using Random
using DataStructures

@enum Face North East South West
@enum NAEnv Ancestral Novel

mutable struct Setting
    seed::Int64
    with_cue::Bool
    max_pop::Int64
    num_cell_x::Int64
    num_cell_y::Int64
    num_env::Int64
    num_layers::Int64
    num_dev::Int64
    num_components::Vector{Int64}
    topology::Vector{OrderedDict{Int64,Float64}}
    omega::Vector{Float64}
    afuncs::Vector{Function}
    state_memory::Vector{Float64}
    env_noise::Float64
    mut_rate::Float64
    conv_dev::Float64
    denv::Float64
    selstrength::Float64
end

const sixpii = 6.0/π
const sqrt3 = √3.0
function lcatan(x)
    sixpii*atan(x/sqrt3)
end

function relu(x)
    max(0, x)
end

default_setting = Setting(13579 # seed
                          , true # with_cue
                          , 200 # max_pop
                          , 1   # num_cell_x
                          , 1   # num_cell_y
                          , 50  # num_env / face
                          , 5   # num_layers
                          , 200 # num_dev 
                          , [200, 200, 200, 200, 200]      # num_components
                          , [OrderedDict(),
                             OrderedDict(1 => 0.02, 2 => 0.02),
                             OrderedDict(2 => 0.02, 3 => 0.02),
                             OrderedDict(3 => 0.02, 4 => 0.02),
                             OrderedDict(4 => 0.02)] # topology
                          , fill(1.0, 5) # omega (use "set_omegas")
                          , [identity, tanh, tanh, tanh, tanh] # afuncs
                          , [0.0, 0.8, 0.0, 0.0, 0.0] # state_memory
                          , 0.05 # env_noise
                          , 0.005 # mut_rate
                          , 1e-5 # conv_dev
                          , 0.5  # magnitude of environmental changes [0,1]
                          , 20.0 # selstrength
                          )

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
