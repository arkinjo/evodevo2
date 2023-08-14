using ArgParse
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
            gps =
                ThreadsX.map(indivs) do indiv
                    g = (genotype(indiv) - genofst) ⋅ dgeno
                    p = (selected_phenotype(indiv, s) - sel0) ⋅ denvs
                    (g,p)
                end
            df = DataFrame(gps)
            rename!(df, :1 => :g, :2 => :p)
            df
        end

        gpdata = DataFrame(gen=Int64[],
                           mg1=Float64[], mp1=Float64[],
                           dg1=Float64[], dp1=Float64[],
                           mg0=Float64[], mp0=Float64[],
                           dg0=Float64[], dp0=Float64[])

        for igen = 1:ngen
            @printf(stderr, "Generation %d\n", igen); flush(stderr)
            pop0 = file[make_pop_name(Ancestral, igen)];
            pop1 = file[make_pop_name(Novel, igen)];

            gp1 = project(pop1.indivs)
            gp0 = project(pop0.indivs)
            
            push!(gpdata, (igen,
                           mean(gp1.g),mean(gp1.p),std(gp1.g), std(gp1.p),
                           mean(gp0.g),mean(gp0.p),std(gp0.g), std(gp0.p)))
            
            scatter(size=(600,600), legend_position=:none,
                    xlims=(-0.05, 1.05), ylims=(-0.05, 1.05),
                    xticks=0:0.1:1, yticks=0:0.1:1,
                    plot_title= @sprintf("%s (gen. %.3d)", s.basename, igen))
            scatter!(gp1.g, gp1.p, label="Novel", markershape=:circle)
            scatter!(gp0.g, gp0.p, label="Anceltral", markershape=:diamond)
            xlabel!("Genotype")
            ylabel!("Phenotype")
            oname= @sprintf("%s_%.3d.pdf", basename, igen)
            Plots.pdf(oname)
        end
        cname= @sprintf("%s.csv", basename)
        CSV.write(cname, gpdata)
            
        scatter(size=(600,600), legend_position=:topleft,
                plot_title= s.basename,
                xlims=(-0.05, 1.05), ylims=(-0.05, 1.05),
                xticks=0:0.1:1, yticks=0:0.1:1)
        scatter!(gpdata.mg1, gpdata.mp1, xerror=gpdata.dg1, yerror=gpdata.dp1,
                 label="Novel", markershape=:circle)
                  
        scatter!(gpdata.mg0, gpdata.mp0, xerror=gpdata.dg0, yerror=gpdata.dp0,
                 label="Anceltral", markershape=:diamond)
        oname= @sprintf("%s_summary.pdf", basename)
        Plots.pdf(oname)
    end
end

main()

