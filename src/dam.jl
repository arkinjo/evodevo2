"""
Correlated environmental cues via the Hopfield network.
"""

struct DAM
    seed::Int64
    ndim::Int64
    patterns::Vector{Vector{Float32}}
end

function DAM(seed, ndim, npat)
    Random.seed!(seed)
    patterns = Vector{Vector{Float32}}()
    for k = 1:npat
        v = rand([-1.0f0, 1.0f0], ndim)
        push!(patterns, v)
    end
    DAM(seed, ndim, patterns)
end

function numpat(hop::DAM)
    length(hop.patterns)
end

function addpat(hop::DAM, pat)
    push!(hop.patterns, pat)
end

function energy(hop::DAM, env)
    -mapreduce(ξ -> exp(ξ ⋅ env/3√hop.ndim), +, hop.patterns)
end

function flip1(hop::DAM, env, i, β)
    ene0 = energy(hop, env)
    env[i] *= -1 # flip
    ene1 = energy(hop, env)
    dene = ene1 - ene0
    if dene > 0.0 && rand() > exp(-β*dene) # Metropolis
        env[i] *= -1 # reject
    end
end

function flip(hop::DAM, env, β)
    for i in shuffle(1:hop.ndim)
        flip1(hop, env, i, β)
    end
end

"batch-sampling all elements at once."
function flip_batch(hop::DAM, env, β)
    dots = map(ξ -> ξ ⋅ env, hop.patterns)
    ene0 = -mapreduce(d -> exp(d/3√hop.ndim), +, dots)
    for (i,e) = enumerate(env)
        ene1 = -mapreduce((d,ξ) -> exp((d - 2e*ξ[i])/3√hop.ndim), +,
                          dots, hop.patterns)
        Δene = ene1 - ene0
        if Δene < 0 || rand() < exp(-β*Δene)
            env[i] *= -1
        end
    end
end

function sample(hop::DAM, n)
    ene = map(_ -> energy(hop, rand([-1.0f0, 1.0f0], hop.ndim)), 1:n)
    μ = mean(ene)
    σ = std(ene, mean=μ)
    nene = map(p -> (energy(hop, p)-μ)/σ, hop.patterns)
    length(hop.patterns), μ, σ, mean(nene), std(nene)
end

function testsize(ndim)
    npats = vcat(map(identity, 50:50:1000), map(identity, 1100:100:2000))
    dams = map(n -> DAM(13, ndim, n), npats)
    DataFrame(map(d -> sample(d, 500), dams), [:npat, :mu, :sig, :zmean, :zsig])
end

"for testing"
function checkhop(hop, ipat, beta, flip)
    env0 = ipat > 0 ?  hop.patterns[ipat] : rand([-1.0f0, 1.0f0], 200)
    ene0 = energy(hop, env0)
    devs = Vector()
    envs = Vector()
    for k = 1:500
        env = copy(env0);
        flip(hop, env, beta)
        dene = energy(hop, env) - ene0
        dev = norm((env - env0)/2, 1)
        push!(envs, env)
        push!(devs, dev)
        @printf("ene0= %f, dene= %f, prob= %f, diff= %f\n",
                ene0, dene, exp(-beta*dene), dev)
    end
    devs, envs
end
