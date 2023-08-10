using JLD2,CodecZlib
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

        "--ngen"
        help = "number of generations per epoch"
        arg_type = Int64
        default = 200
        
        "restart"
        help = "restart JDL2 file"
        required = true
    end

    return parse_args(s)
end


function main()
    parsed_args = parse_commandline()
    nepoch = parsed_args["nepoch"]
    ngen = parsed_args["ngen"]
    seed = parsed_args["seed"]
    resfile = parsed_args["restart"]
    s, envs, pop =
        jldopen(resfile, "r") do file
            s = file["setting"]
            envs = file["envs"]
            pop = file["pop0"]
            s,envs, pop
        end
    open(s.basename * "_test.dat", "w") do log
        println(log, "#basename= $s")
        flush(log)
        envs0 = copy(envs)
        for iepoch = 1:nepoch
            trajfile = @sprintf("%s_traj%.2d.jld2", s.basename, iepoch)
            s.seed += iepoch
            envs1 = change_envS(envs0, s)
            jldopen(trajfile, "w") do traj
                traj["setting"] = s
                traj["epoch"] = iepoch
                traj["envs0"] = envs0
                traj["envs1"] = envs1
                pop = evolve(TestMode, iepoch, ngen, pop, envs0, envs1, log,
                             traj, s)
            end
            envs0 = envs1
            flush(log)
        end
    end
end

main()

