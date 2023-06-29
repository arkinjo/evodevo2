using LinearAlgebra
using Statistics

struct Individual
    id::Int64
    mom_id::Int64
    dad_id::Int64
    genome::Genome
    cells::Dict{Tuple{Int64,Int64},Union{Cell,Environment}}
    ndev::Vector{Float64}
    mismatch::Vector{Float64}
    fitness::Vector{Float64}
end


function Individual(id, dad, mom, genome::Genome, s::Setting)
    cells = Dict{Tuple{Int64,Int64},Union{Cell,Environment}}()
    for i = 1:s.num_cell_x
        for j = 1:s.num_cell_y
            cells[i,j] = Cell(s.num_components)
        end
    end
    Individual(id, dad, mom, genome, cells, [0], [0], [0, 0])
end

function Individual(id::Int64, s::Setting)
    genome = Genome(s)
    Individual(id, 0, 0, genome, s)
end

function phenotype(indiv::Individual, s::Setting)
    pheno = Vector()
    for i = 1:s.num_cell_x, j = 1:s.num_cell_y
        p = phenotype(indiv.cells[i,j])
        pheno = vcat(pheno, p)
    end
    pheno
end

function get_selected_phenotype(indiv::Individual, s::Setting)
    f = Vector{Float64}()
    for i = 1:s.num_cell_x
        f = vcat(f, get_face(indiv.cells[i, s.num_cell_y], North))
    end
    f
end
    
function set_fitness(indiv::Individual, envs::EnvironmentS, s::Setting)
    e = get_selecting_envs(envs, s)
    f = get_selected_phenotype(indiv, s)
    tdev = (e-f) .|> abs |> mean
    indiv.mismatch[] = tdev
    if indiv.ndev[] < s.num_dev
        wdev = max(0, indiv.ndev[] - 100)/20.0
        indiv.fitness[1] = exp(-tdev*s.selstrength - wdev)
    else
        indiv.fitness[1] = 0.0
    end
end

function set_envs(indiv::Individual, envs::EnvironmentS, s::Setting)
    for (k, env) in envs
        indiv.cells[k] = get_cue(env, s)
    end
end

function develop(indiv::Individual, envs::EnvironmentS, s::Setting) ::Bool
    ndev = 0
    while ndev < s.num_dev
        ndev += 1
        for i = 1:s.num_cell_x, j = 1:s.num_cell_y
            north = get_face(indiv.cells[i,j+1], South)
            east = get_face(indiv.cells[i+1,j], West)
            south = get_face(indiv.cells[i,j-1], North)
            west = get_face(indiv.cells[i-1,j], East)
            indiv.cells[i,j].states[1] = vcat(north, east, south, west)
            DevStep(ndev, indiv.cells[i,j], indiv.genome, s)
        end
        dev = get_dev(indiv, s)
        if ndev > 1 && dev < s.conv_dev
            break
        end
    end
    set_fitness(indiv, envs, s)
    indiv.ndev[] = ndev
    (ndev < s.num_dev)
end

function get_dev(indiv::Individual, s::Setting)
    pvar = 0.0
    n = 0.0
    for i = 1:s.num_cell_x, j = 1:s.num_cell_y
        pvar += sum(indiv.cells[i,j].pvar)
        n += length(indiv.cells[i,j].pvar)
    end
    pvar/n
end

function mate(ind1::Individual, ind2::Individual)
    (geno1, geno2) = mate(ind1.genome, ind2.genome)
end
