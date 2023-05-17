struct Body
    genome::Genome
    cells::Vector{Cell}
end    

struct Individual
    id::Int64
    mom_id::Int64
    dad_id::Int64
    bodies::Vector{Body}
