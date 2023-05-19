using Random

Environment = Vector{Float64}

function make_environment(seed, n)
    Random.seed!(seed)
    rand([-1, 1], n)
end

function get_face(env, i)
    n = length(env) ÷ 4
    env[(i-1)*n+1:i*n]
end
