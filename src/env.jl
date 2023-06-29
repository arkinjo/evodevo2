using Random

struct Environment
    p::Vector{Float64}
end

"""
Environments is the surrounding environment of individuals.
"""
EnvironmentS = Dict{Tuple{Int64, Int64}, Environment}

"""
    make_environments(s::Setting)
"""
function make_environments(s::Setting) ::EnvironmentS
    Random.seed!(s.seed)
    envs = Dict{Tuple{Int64,Int64},Environment}()
    for i = 1:s.num_cell_x
        envs[i, 0] = Environment(rand([-1.0, 1.0], s.num_env))
        envs[i, s.num_cell_y+1] = Environment(rand([-1.0, 1.0], s.num_env))
    end
    for j = 1:s.num_cell_y
        envs[0, j] = Environment(rand([-1.0, 1.0], s.num_env))
        envs[s.num_cell_x + 1, j] = Environment(rand([-1.0, 1.0], s.num_env))
    end
    envs
end

function get_face(env::Environment, _)
    env.p
end

function get_cue(env::Environment, s::Setting)
    cue = s.with_cue ? copy(env.p) : fill(1.0, s.num_env)
    m = Int64(round(s.env_noise*length(cue)))
    r = randperm(length(cue))
    for i in r[1:m]
#        cue[i] = s.with_cue ? -cue[i] : rand([-1.0, 1.0])
        cue[i] = -cue[i]
    end
    Environment(cue)
end

function change_env(env::Environment, s::Setting)
    m = Int(round(length(env.p)*s.denv))
    r = randperm(s.num_env)
    p = copy(env.p)
    for i = r[1:m]
        p[i] *= -1
    end
    Environment(p)
end

function change_envS(envs::EnvironmentS, seed, s::Setting)
    Random.seed!(seed)
    nenvs = Dict{Tuple{Int64,Int64},Environment}()
    a = 0
    for (k,v) in envs
        a += 1
        nenvs[k] = change_env(v, s)
    end
    nenvs
end

function get_selecting_envs(envs::EnvironmentS, s::Setting)
    e = Vector{Float64}()
    for i = 1:s.num_cell_x
        e = vcat(e, get_face(envs[i, s.num_cell_y+1], South))
    end
    e
end
