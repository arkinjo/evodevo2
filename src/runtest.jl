using JLD2,CodecBzip2
using ArgParse

include("EvoDevo2.jl")

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--seed"
        help = "random seed (13579)"
        arg_type = Int64
        default = 13579

        "--outdir"
        help = "output directory for trajectory JLD2 files"
        arg_type = String
        default = "."

        "--nepoch"
        help = "number of epochs"
        arg_type = Int64
        default = 1

        "--ngen"
        help = "number of generations per epoch"
        arg_type = Int64
        default = 200

        "--denv"
        help = "magnitude of environmental change [0.0 ~ 0.5]"
        arg_type = Float64
        default = 0.5

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
    outdir = parsed_args["outdir"]
    s, envs, pop =
        jldopen(resfile, "r") do file
            s = file["setting"]
            envs = file["envs"]
            pop = file["pop"]
            s,envs,pop
        end
    s.denv = parsed_args["denv"]

    open(s.basename * "_test.dat", "w") do log
        println(log, "#basename= $s")
        flush(log)
        envs0 = copy(envs)
        for iepoch = 1:nepoch
            trajfile = @sprintf("%s/%s_traj%.2d.jld2",
                                outdir, s.basename, iepoch)
            s.seed += iepoch
            envs1 = change_envS(envs0, s)
            jldopen(trajfile, "w"; compress=Bzip2Compressor()) do traj
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

