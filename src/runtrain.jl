using JLD2
using ArgParse
using Printf

include("EvoDevo2.jl")

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--seed"
        help = "random seed (13579)"
        arg_type = Int64
        default = 13579

        "--outdir"
        help = "output directory for restart JLD2 files"
        arg_type = String
        default = "."

        "--nepoch"
        help = "number of epochs"
        arg_type = Int64
        default = 1

        "--max_pop"
        help = "population size"
        arg_type = Int64
        default = 200
        
    end

    return parse_args(s)
end


function set_full(seed)
    s = default_setting("Full", 5, seed)
    set_omegas(s)
end

function set_nohier(seed)
    s = default_setting("NoHier", 3, seed)
    s.num_components[2] = 600
    s.topology[2][1] *= (1.0/3.0)
    s.topology[2][2] *= (5.0/9.0)
    s.topology[3][2] *= (1.0/3.0)
    set_omegas(s)
end

function set_nocue(seed)
    s = default_setting("NoCue", 5, seed)
    s.with_cue = false
    set_omegas(s)
end

function set_nodev(seed)
    s = default_setting("NoDev", 5, seed)
    s.num_dev = 1
    set_omegas(s)
end
    
function main()
    parsed_args = parse_commandline()
    max_pop = parsed_args["max_pop"]
    nepoch = parsed_args["nepoch"]
    outdir = parsed_args["outdir"]
    seed = parsed_args["seed"]

    models = ["Full" => set_full(seed),
              "NoHier" => set_nohier(seed),
              "NoCue" => set_nocue(seed),
              "NoDev" => set_nodev(seed)
              ]

    for (model, s) in models
        @printf(stderr, "Model: %s\n", model)
        basename = s.basename
        s.basename = basename
        s.max_pop = max_pop
        logfile = @sprintf("%s/%s.dat", outdir, basename)
        open(logfile * ".dat", "w") do log
            println(log,"epoch\tgen\tmis1\tfit1\tndev1\tali1\tpar1")
            flush(log)
            envs, pop = train_epochs(nepoch, 200, log, s)
            ofilename = @sprintf("%s/%s_train.jld2", outdir, basename)
            jldopen(ofilename, "w") do file
                file["setting"] = s
                file["envs"] = envs
                file["pop"] = pop
            end
        end
    end
end

main()
