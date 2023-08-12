using JLD2,CodecBzip2
using ArgParse
using Statistics
using ThreadsX
using Plots
using Printf

include("EvoDevo2.jl")

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--outdir"
        help = "output directory for G-P plot"
        arg_type = String
        default = "gp"

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
        sel0 = selecting_envs(envs0,s)
        sel1 = selecting_envs(envs1,s)
        denvs = sel1 - sel0
        denvs /= dot(denvs, denvs)

        popfst = file[make_pop_name(Ancestral,1)];
        poplst = file[make_pop_name(Novel, ngen)];

        genofst = Statistics.mean(get_geno_vecs(popfst); dims=2)
        genolst = Statistics.mean(get_geno_vecs(poplst); dims=2)
        dgeno = genolst - genofst
        dgeno /= dot(dgeno, dgeno)
        
        function project(indivs)
            ThreadsX.map(indivs) do indiv
                g = (genotype(indiv) - genofst) ⋅ dgeno
                p = (selected_phenotype(indiv, s) - sel0) ⋅ denvs
                [g,p]
            end
        end
        mgs0 = Vector(); dgs0 = Vector()
        mgs1 = Vector(); dgs1 = Vector()
        mps0 = Vector(); dps0 = Vector()
        mps1 = Vector(); dps1 = Vector()

        for igen = 1:ngen
            pop0 = file[make_pop_name(Ancestral, igen)]
            pop1 = file[make_pop_name(Novel, igen)]

            gp0 = project(pop0.indivs)
            (mg0,mp0) = mean(gp0)
            (dg0,dp0) = std(gp0)
            push!(mgs0,mg0); push!(dgs0, dg0)
            push!(mps0,mp0); push!(dps0, dp0)
            gp1 = project(pop1.indivs)
            (mg1,mp1) = mean(gp1)
            (dg1,dp1) = std(gp1)
            push!(mgs1,mg1); push!(dgs1, dg1)
            push!(mps1,mp1); push!(dps1, dp1)

            @printf("%d\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",
                    igen, mg1, mp1, dg1, dp1, mg0, mp0, dg0, dp0)

            g0 = map(x->x[1], gp0); p0 = map(x->x[2], gp0)
            g1 = map(x->x[1], gp1); p1 = map(x->x[2], gp1)
            scatter(size=(600,600), legend_position=:none,
                    xlims=(-0.05, 1.05), ylims=(-0.05, 1.05),
                    xticks=0:0.1:1, yticks=0:0.1:1,
                    plot_title= @sprintf("%s (gen. %.3d)", s.basename, igen))
            scatter!(g1,p1, label="Novel", markershape=:circle)
            scatter!(g0,p0, label="Anceltral", markershape=:diamond)
            xlabel!("Genotype")
            ylabel!("Phenotype")
            oname= @sprintf("%s_%.3d.pdf", basename, igen)
            Plots.pdf(oname)
        end
        scatter(size=(600,600), legend_position=:topleft,
                plot_title= s.basename,
                xlims=(-0.05, 1.05), ylims=(-0.05, 1.05),
                xticks=0:0.1:1, yticks=0:0.1:1)
        scatter!(mgs1,mps1, xerror=dgs1, yerror=dps1,
                 label="Novel", markershape=:circle)
        scatter!(mgs0,mps0, xerror=dgs0, yerror=dps0,
                 label="Anceltral", markershape=:diamond)
        oname= @sprintf("%s_summary.pdf", basename)
        Plots.pdf(oname)
    end
end

main()

