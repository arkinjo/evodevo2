using Random
using Distributions
using SparseArrays
using DataStructures
using JLD2
using Printf

"Each genome matrix is a sparse matrix"
GenomeMat = SparseMatrixCSC{Float32,Int64}

"""
F is an array of feedforward matrices.
B is an array of dictionaries of feedback matrices
"""
struct Genome
    B::Vector{OrderedDict{Int64,GenomeMat}}
end

struct Mutation
    mats::Vector{Any}
    cats::Categorical{Float32, Vector{Float32}}
    num_mut::Poisson{Float32}
end

function Mutation(s::Setting)
    mats = Vector()
    probs = Vector{Float32}()
    tot = 0.0
    for (l, m) in enumerate(s.topology)
        for (k, d) in m
            push!(mats, (l,k))
            nc = s.num_components[l]*s.num_components[k]
            tot += nc
            push!(probs, nc)
        end
    end
    for k in keys(probs)
        probs[k] /= tot
    end
    cats = Categorical(probs)
    num_mut = Poisson(s.mut_rate*tot)
    Mutation(mats, cats, num_mut)
end

"""
    make_genome_mat(m, n, density)
Make a single genome matrix, `m`-by-`n` sparse matrix with `density`.
"""
function make_genome_mat(m, n, density)
    I = Vector{Int64}()
    J = Vector{Int64}()
    V = Vector{Float32}()
    for i = 1:m, j = 1:n
        if rand() < density
            push!(I,i)
            push!(J,j)
            push!(V,rand([-1.0, 1.0]))
        end
    end
    sparse(I,J,V, m, n)
end

"""
Genome(s::Setting)
Make the genome of an individual with `nlayers` layers with `ncomps` components and the feedback loops given by `fbloops`.
"""
function Genome(s::Setting)
    b = Array{OrderedDict{Int64,GenomeMat}}(undef, s.num_layers)
    for l = 1:s.num_layers
        d = OrderedDict{Int64,GenomeMat}()
        for (k, ρ) in s.topology[l]
            d[k] = make_genome_mat(s.num_components[l],s.num_components[k], ρ)
        end
        b[l] = d
    end
    Genome(b)
end

function mutate(A::GenomeMat, density::Float32, mut_rate::Float32)
    (m,n) = size(A)
    poisson = Poisson(mut_rate*m*n)
    nmut = rand(poisson)
    d2 = density*0.5
    for k = 1:nmut
        i = rand(1:m)
        j = rand(1:n)
        q = rand()
        A[i,j] = 0.0
        if q < d2
            A[i,j] = -1.0
        elseif q < density
            A[i,j] = 1.0
        end
    end
    dropzeros!(A)
end

function mutate1(A::GenomeMat, density::Float32)
    (m,n) = size(A)
    i = rand(1:m)
    j = rand(1:n)

    if rand() >= density
        A[i,j] = 0.0
    else
        A[i,j] = rand((-1.0, 1.0))
    end
    A
end

function mutate(genome::Genome, muts::Mutation, s::Setting)
    nmut = rand(muts.num_mut)
    for _ = 1:nmut
        (l,k) = muts.mats[rand(muts.cats)]
        mutate1(genome.B[l][k], s.topology[l][k])
    end
    for l = 2:s.num_layers
        for (k,b) in genome.B[l]
            dropzeros!(genome.B[l][k])
        end
    end
end

function mate(mat1::GenomeMat, mat2::GenomeMat)
    (m,n) = size(mat1)
    nmat1 = spzeros(m,n)
    nmat2 = spzeros(m,n)
    for i = 1:m
        r = rand()
        if r < 0.5
            nmat1[i,:] = copy(mat1[i,:])
            nmat2[i,:] = copy(mat2[i,:])
        else
            nmat2[i,:] = copy(mat1[i,:])
            nmat1[i,:] = copy(mat2[i,:])
        end
    end
    dropzeros!(nmat1)
    dropzeros!(nmat2)
    nmat1, nmat2
end

function mate(geno1::Genome, geno2::Genome)
    B1 = Vector{OrderedDict{Int64,GenomeMat}}()
    B2 = Vector{OrderedDict{Int64,GenomeMat}}()
    for l in 1:length(geno1.B)
        Bk1 = OrderedDict{Int64, GenomeMat}()
        Bk2 = OrderedDict{Int64, GenomeMat}()
        for (k, b1) in geno1.B[l]
            b2 = geno2.B[l][k]
            nb1, nb2 = mate(b1, b2)
            Bk1[k] = nb1
            Bk2[k] = nb2
        end
        push!(B1, Bk1)
        push!(B2, Bk2)
    end
    Genome(B1), Genome(B2)
end

function vectorize(geno::Genome)
    len = 0
    for d in geno.B
        for M in values(d)
            m,n = size(M)
            len += m*n
        end
    end
    vec = zeros(len)
    for d in geno.B
        for M in values(d)
            i = 1
            for x in values(M)
                vec[i] = x
                i += 1
            end
        end
    end
    vec
end

