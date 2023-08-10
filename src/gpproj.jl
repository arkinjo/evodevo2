using JLD2,CodecZlib
using ArgParse
using Statistics

include("EvoDevo2.jl")

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--outdir"
        help = "output directory for G-P plot"
        arg_type = String
        default = "."

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
    dir= parsed_args["outdir"]
    trajfile = parsed_args["traj"]
    jldopen(trajfile, "r") do file
        s = file["setting"]
        epoch = file["epoch"]
        basename = @sprintf("%s/%s_gp%.2d", dir, s.basename, epoch)
        envs0 = file["envs0"]
        envs1 = file["envs1"]
        popfst = file["pop0_001"];
        lst = @sprintf("pop1_%.3d", ngen)
        poplst = file[lst];

        genofst = Statistics.mean(get_geno_vecs(popfst); dims=2)
        genolst = Statistics.mean(get_geno_vecs(poplst); dims=2)
        dgeno = genolst - genofst
        dgeno /= norm(dgeno)^2

        sel0 = selecting_envs(envs0,s)
        sel1 = selecting_envs(envs1,s)
        denvs = sel1 - sel0
        denvs /= norm(denvs)^2
        
        function project(indivs)
            map(indivs) do indiv
                g = (genotype(indiv) - genofst) ⋅ dgeno
                p = (selected_phenotype(indiv, s) - sel0) ⋅ denvs
                [g,p]
            end
        end
        for igen = 1:ngen
            name0 = @sprintf("pop0_%.3d", igen)
            name1 = @sprintf("pop1_%.3d", igen)
            pop0 = file[name0]
            pop1 = file[name1]
            gp0 = project(pop0.indivs)
            gp1 = project(pop1.indivs)
            (mg0,mp0) = mean(gp0)
            (dg0,dp0) = std(gp0)
            
            (mg1,mp1) = mean(gp1)
            (dg1,dp1) = std(gp1)

            @printf("%d\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",
                    igen, mg1, mp1, dg1, dp1, mg0, mp0, dg0, dp0)

            oname= @sprintf("%s_%.3d.dat", basename, igen)
            open(oname, "w") do fh
                for (i,((g1,p1), (g0, p0))) in enumerate(zip(gp1, gp0))
                    @printf(fh, "%d\t%d\t%e\t%e\t%e\t%e\n",
                            igen, i, g1, p1, g0, p0)
                end
            end
        end
    end
end

main()

