struct Cell
    id::Int64
    states::Vector{Vector{Float64}}
    pave::Vector{Float64}
    pvar::Vector{Float64}
    perr::Float64
end

function Cell(id, nlayers, ncomps)
    pave = zeros(ncomps[end])
    pvar = zeros(ncomps[end])
    states = Vector{Array{Float64}}()
    for n in ncomps
        push!(states, ones(n))
    end
    Cell(id, states, pave, pvar, Inf)
end

phenotype(cell) = cell.states[end]

struct Body
    cells
end

