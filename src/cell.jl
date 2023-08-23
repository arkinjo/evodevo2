using Statistics

struct Cell
    states::Vector{Vector{Float32}}
    pave::Vector{Float32}
    pvar::Vector{Float32}
end

function Cell(ncomps::Vector{Int64})
    states = map(n -> ones(Float32, n), ncomps)
    pave = zeros(length(states[end]))
    pvar = zeros(length(states[end]))
    Cell(states, pave, pvar)
end

function (cell::Cell)(k)
    cell.states[k]
end

phenotype(cell::Cell) = cell.pave # cell.states[end]

function get_face(cell::Cell, face::Face)
    p = phenotype(cell)
    n = length(p) ÷ 4
    i = Int64(face)
    p[i*n+1:(i+1)*n]
end

function update_pave(nstep, cell::Cell, s::Setting)
    mstep = (nstep < 5) ? Float32(nstep) : 5.0
    α = 2.0/(mstep + 1.0)

    for (i,(p,a,v)) in enumerate(zip(cell.states[end], cell.pave, cell.pvar))
        d = p - a
        incr = α * d
        cell.pave[i] += incr
        cell.pvar[i] = (1 - α)*(v + d*incr)
    end
    mean(cell.pvar)
end

function DevStep(istep, cell, genome, s::Setting)
    for l = 2:s.num_layers
        ϕ = s.state_memory[l]*cell.states[l]
        for (k,mat) in genome.B[l]
            if l == 2 && k == 1
                ϕ += mat * (cell.states[k] - phenotype(cell))
            else
                ϕ += mat * cell.states[k]
            end
        end
        cell.states[l] = s.afuncs[l].(ϕ * s.omega[l]) 
    end
    update_pave(istep, cell, s)
end

function check_vecs(v1, v2)
    for (i, (a,b)) in enumerate(zip(v1, v2))
        if abs(a-b) > 0
            println("$i : $a $b")
        end
    end
end
