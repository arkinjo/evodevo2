using JLD2,CodecZlib
using ArgParse
using Statistics

include("EvoDevo2.jl")

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--out"
        help = "output restart JLD2 file"
        arg_type = String
        default = "proj.dat"

        "--ngen"
        help = "number of generations"
        arg_type = Int64
        default = 200
        
        "traj"
        help = "trajectory JLD2 file"
        arg_type = String
        required = true
        
    end

    return parse_args(s)
end


function main()
    parsed_args = parse_commandline()
    ngen = parsed_args["ngen"]
    trajfile = parsed_args["traj"]
    jldopen(trajfile, "r") do file
        s = file["setting"]
        envs0 = file["envs0"]
        envs1 = file["envs1"]
        popfst = file["pop0_001"];
        lst = @sprintf("pop1_%.3d", ngen)
        poplst = file[lst];

        genofst = Statistics.mean(get_geno_vecs(popfst); dims=2)
        genolst = Statistics.mean(get_geno_vecs(poplst); dims=2)
        dgeno = genolst - genofst
        dgeno /= norm(dgeno)*norm(dgeno)
        
        sel0 = get_selecting_envs(envs0,s)
        sel1 = get_selecting_envs(envs1,s)
        denvs = sel1 - sel0
        denvs /= norm(denvs)*norm(denvs)
        for i = 1:ngen
            name0 = @sprintf("pop0_%.3d", ngen)
            name1 = @sprintf("pop1_%.3d", ngen)
            pop0 = file[name0]
            pop1 = file[name1]
        end
    end
end

main()

