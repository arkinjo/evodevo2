using JLD2
using ArgParse

include("EvoDevo2.jl")

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--seed"
        help = "random seed (13579)"
        arg_type = Int64
        default = 13579
        
        "--nepoch"
        help = "number of epochs"
        arg_type = Int64
        default = 1

        "--max_pop"
        help = "population size"
        arg_type = Int64
        default = 200
        
        "ext"
        help = "basename extension for output files"
        default = "00"
        required = false
    end

    return parse_args(s)
end


function set_full(seed)
    s = default_setting("full", 5, seed)
    set_omegas(s)
end

function set_nohier(seed)
    s = default_setting("nohier", 3, seed)
    s.num_components[2] = 600
    set_omegas(s)
end

function set_nocue(seed)
    s = default_setting("nocue", 5, seed)
    s.with_cue = false
    set_omegas(s)
end

function set_nodev(seed)
    s = default_setting("nodev", 5, seed)
    s.num_dev = 1
    set_omegas(s)
end
    
function main()
    parsed_args = parse_commandline()
    max_pop = parsed_args["max_pop"]
    nepoch = parsed_args["nepoch"]
    seed = parsed_args["seed"]
    ext = parsed_args["ext"]

    models = ["full" => set_full(seed),
              "nohier" => set_nohier(seed),
              "nocue" => set_nocue(seed),
              "nodev" => set_nodev(seed)
              ]

    Threads.@threads for (model, s) in models
        basename = s.basename * ext
        s.max_pop = max_pop
        open(basename * ".dat", "w") do log
            println(log, "#basename= $s")
            flush(log)
            envs, pop = train_epochs(nepoch, 200, log, s)
            jldopen(basename * "_train.jld2", "w") do file
                file["setting"] = s
                file["envs"] = envs
                file["pop0"] = pop
            end
        end
    end
end

main()
