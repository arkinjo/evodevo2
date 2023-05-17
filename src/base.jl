using Random

struct Setting
    num_layers::Int64
    num_components::Vector{Int64}
    feedback_loops::Vector{Vector{Int64}}
end

default_setting = Setting(5, [10, 10, 10, 10, 10],
                          [[], [3, 5], [4], [], []])

