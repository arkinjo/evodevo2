using Random

struct Setting
    num_cell_x::Int64
    num_cell_y::Int64
    num_env::Int64
    num_layers::Int64
    num_components::Vector{Int64}
    feedback_loops::Vector{Vector{Int64}}
end

default_setting = Setting(3, 3, 5, 5, [20, 20, 20, 20, 4*num_env],
                          [[], [3, 5], [4], [], []])

