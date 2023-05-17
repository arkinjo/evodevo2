using Random
using SparseArrays
using JLD2

"Each genome matrix is a sparse matrix"
GenomeMat = SparseMatrixCSC{Float64,Int64}

"""
F is an array of feedforward matrices.
B is an array of dictionaries of feedback matrices
"""
struct Genome
    F::Vector{GenomeMat}
    B::Vector{Dict{Int64,GenomeMat}}
end

"""
    make_genome_mat(m, n, density)
Make a single genome matrix, `m`-by-`n` sparse matrix with `density`.
"""
function make_genome_mat(m, n, density)
    I = Vector{Int64}()
    J = Vector{Int64}()
    V = Vector{Float64}()
    for i = 1:m, j = 1:n
        if rand() < density
            push!(I,i)
            push!(J,j)
            push!(V,2.0 * Random.bitrand()[1] - 1.0)
        end
    end
    sparse(I,J,V)
end

"""
    make_genome(nlayers, ncomps, fbloops)
Make the genome of an individual with `nlayers` layers with `ncomps` components and the feedback loops given by `fbloops`.
"""
function make_genome(nlayers, ncomps, fbloops)
    f = Vector{GenomeMat}()
    push!(f, spzeros(0, ncomps[1]))
    for i = 2:nlayers
        push!(f,make_genome_mat(ncomps[i-1],ncomps[i], 0.1))
    end
    b = Array{Dict{Int64,GenomeMat}}(undef, nlayers)
    for i = 2:nlayers
        d = Dict{Int64,GenomeMat}()
        for j in fbloops[i]
            d[j] = make_genome_mat(ncomps[i],ncomps[j], 0.1)
        end
        b[i] = d
    end
    Genome(f, b)
end
