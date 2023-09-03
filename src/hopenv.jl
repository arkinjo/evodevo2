"""
Correlated environmental cues via the Hopfield network.
"""
using Random
using LinearAlgebra
using Printf

struct HopEnv
    seed::Int64
    ndim::Int64
    npat::Int64
    patterns::Vector{Vector{Float32}}
    W::Matrix{Float32}
end

function HopEnv(seed, ndim, npat)
    Random.seed!(seed)
    patterns = Vector{Vector{Float32}}()
    W = zeros(Float32, ndim, ndim)
    for k = 1:npat
        v = rand([-1.0f0, 1.0f0], ndim)
        push!(patterns, v)
        W += v*v' 
    end
    W[diagind(W)] .= 0.0
    W /= npat
    HopEnv(seed, ndim, npat, patterns, W)
end

function flip1(hop::HopEnv, env, i, β)
    w = hop.W[i,:] ⋅ env
    if rand() < exp(-β*env[i]*w) # Metropolis
            env[i] *= -1
    end
end

function flipn(hop::HopEnv, n, β, env)
    for i in shuffle(1:n)
        flip1(hop, env, i, β)
    end
end

function DAMenergy(hop::HopEnv, env)
    energy = mapreduce(ξ -> exp(ξ ⋅ env/2√hop.ndim), +, hop.patterns)
    energy *= -1/hop.npat
end

function DAMflip1(hop::HopEnv, env, i, β)
    ene0 = DAMenergy(hop, env)
    env[i] *= -1
    ene1 = DAMenergy(hop, env)
    dene = ene1 - ene0
    if dene > 0.0 && rand() > exp(-β*dene) # Metropolis
        env[i] *= -1 # reject
    end
end

function DAMflip(hop::HopEnv, env, β)
    for i in shuffle(1:hop.ndim)
        DAMflip1(hop, env, i, β)
    end
end

"batch-sampling all elements at once."
function DAMflip_batch(hop::HopEnv, env, β)
    dots = map(ξ -> ξ ⋅ env, hop.patterns)
    ene0 = -mapreduce(d -> exp(d/2√hop.ndim), +, dots)/hop.npat
    for (i,e) = enumerate(env)
        ene1 = -mapreduce((d,ξ) -> exp((d - 2e*ξ[i])/2√hop.ndim), +,
                          dots, hop.patterns)/hop.npat
        Δene = ene1 - ene0
        if Δene < 0 || rand() < exp(-β*Δene)
            env[i] *= -1
        end
    end
end

"for testing"
function checkhop(hop, ipat, beta, flip)
    env0 = ipat > 0 ?  hop.patterns[ipat] : rand([-1.0f0, 1.0f0], 200)
    ene0 = DAMenergy(hop, env0)
    devs = Vector()
    envs = Vector()
    for k = 1:500
        env = copy(env0);
        flip(hop, env, beta)
        dene = DAMenergy(hop, env) - ene0
        dev = norm((env - env0)/2, 1)
        push!(envs, env)
        push!(devs, dev)
        @printf("ene0= %f, dene= %f, prob= %f, diff= %f\n",
                ene0, dene, exp(-beta*dene), dev)
    end
    devs, envs
end
