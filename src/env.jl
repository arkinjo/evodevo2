using Random

struct Environment
    p::Vector{Float32}
end

"""
Environments is the surrounding environment of individuals.
"""
struct EnvironmentS
    p::Vector{Float32}
    dam::DAM
end

"""
    make_environments(s::Setting)
"""
function make_environments(s::Setting) ::EnvironmentS
    ndim = s.num_env * (s.num_cell_x + s.num_cell_y) * 2
    dam = DAM(s.seed, ndim, 100)
    p = copy(dam.patterns[1])
    EnvironmentS(p, dam)
end

function get_face(env::Environment, _)
    env.p
end

function get_cues(envs::EnvironmentS, s::Setting)
    cue = s.with_cue ? copy(envs.p) : ones(Float32, length(envs.p))
    m = Int64(round(s.env_noise*length(cue)))
    r = randperm(length(cue))
    for i in r[1:m]
        cue[i] = -cue[i]
    end
    envs_of_vec(cue, s)
end

function change_envS(envs::EnvironmentS, s::Setting)
    Random.seed!(s.seed)
    m = Int(round(length(envs.p)*s.denv))
    r = randperm(length(envs.p))
    p = copy(envs.p)
    for i = r[1:m]
        p[i] *= -1
    end
    EnvironmentS(p, envs.dam)
end

function selecting_envs(envs::EnvironmentS, s::Setting)
    dict = envs_of_vec(envs.p, s)
    e = Vector{Float32}()
    for i = 1:s.num_cell_x
        e = vcat(e, get_face(dict[i, s.num_cell_y+1], South))
    end
    e
end

function flatten(envs::EnvironmentS, s::Setting)
    v = Vector{Float32}()
    for i = 1:s.num_cell_x
        v = vcat(v,envs[i, 0].p)
        v = vcat(v, envs[i, s.num_cell_y+1].p)
    end
    for j = 1:s.num_cell_y
        v = vcat(v, envs[0, j].p)
        v = vcat(v, envs[s.num_cell_x + 1, j].p)
    end
    v
end

function envs_of_vec(v, s::Setting)
    envs = Dict{Tuple{Int64,Int64},Environment}()
    k = 0
    for i = 1:s.num_cell_x
        envs[i, 0] = Environment(v[k*s.num_env+1:(k+1)*s.num_env])
        k += 1
        envs[i, s.num_cell_y+1] = Environment(v[k*s.num_env+1:(k+1)*s.num_env])
        k += 1 
    end
    for j = 1:s.num_cell_y
        envs[0, j] = Environment(v[k*s.num_env+1:(k+1)*s.num_env])
        k += 1
        envs[s.num_cell_x + 1, j] = Environment(v[k*s.num_env+1:(k+1)*s.num_env])
        k += 1
    end
    envs
end

