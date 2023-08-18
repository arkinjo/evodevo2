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

        "--out"
        help = "output CSV file"
        arg_type = String
        default = "gp/gp.csv"

        "--ngen"
        help = "number of generations"
        arg_type = Int64
        default = 200

        "traj"
        nargs = '+'
        help = "trajectory JLD2 files"
        arg_type = String
    end

    return parse_args(s)
end

function mkdf()
    DataFrame(basename=String[], epoch=Int64[], gen=Int64[],
              mg1=Float64[], mp1=Float64[],
              dg1=Float64[], dp1=Float64[],
              gvar1=Float64[], pvar1=Float64[], mis1in1=Float64[], mis1in0=Float64[],
              mg0=Float64[], mp0=Float64[],
              dg0=Float64[], dp0=Float64[],
              gvar0=Float64[], pvar0=Float64[], mis0in1=Float64[], mis0in0=Float64[])
end

function proctraj(dir, ngen, trajfile)
    gpdata = mkdf()
    jldopen(trajfile, "r") do file
        s = file["setting"]
        epoch = file["epoch"]
        @printf(stderr, "%s: Epoch: %d (thread: %d)\n", s.basename, epoch, Threads.threadid()); flush(stderr)
        
        basename = @sprintf("%s/%s_gp%.2d", dir, s.basename, epoch)
        envs0 = file["envs0"]
        envs1 = file["envs1"]
        sel0 = selecting_envs(envs0,s)
        sel1 = selecting_envs(envs1,s)
        denvs = sel1 - sel0
        denvs /= dot(denvs, denvs)

        popfst = file[make_pop_name(Ancestral,1)];
        poplst = file[make_pop_name(Novel, ngen)];

        npop = length(popfst.indivs)
        
        genofst = Statistics.mean(get_geno_vecs(popfst); dims=2)
        genolst = Statistics.mean(get_geno_vecs(poplst); dims=2)
        dgeno = genolst - genofst
        dgeno /= dot(dgeno, dgeno)

        function project(geno,pheno)
            gps =
                ThreadsX.map(1:npop) do i
                    g = (geno[:,i] - genofst) ⋅ dgeno
                    p = (pheno[:,i] - sel0) ⋅ denvs
                    m1 = norm(pheno[:,i] - sel1, 1)
                    m0 = norm(pheno[:,i] - sel0, 1)
                    (g,p,m1,m0)
                end
            df = DataFrame(gps)
            rename!(df, :1 => :g, :2 => :p, :3 => :m1, :4 => :m0)
            df
        end

        for igen = 1:ngen
            @printf(stderr, "\tGeneration %d (thread: %d)\n", igen, Threads.threadid()); flush(stderr)
            pop1 = file[make_pop_name(Novel, igen)];
            pop0 = file[make_pop_name(Ancestral, igen)];
            
            geno1 = get_geno_vecs(pop1)
            pheno1 = get_selected_pheno_vecs(pop1, s)
            
            geno0 = get_geno_vecs(pop0)
            pheno0 = get_selected_pheno_vecs(pop0, s)

            gp1 = project(geno1, pheno1)
            gp0 = project(geno0, pheno0)
            
            sp = scatter(size=(600,600), legend_position=:none,
                    xlims=(-0.05, 1.05), ylims=(-0.05, 1.05),
                    xticks=0:0.1:1, yticks=0:0.1:1,
                    xlabel="Genotype", ylabel="Phenotype",
                    plot_title= @sprintf("%s (gen. %.3d)", s.basename, igen))
            scatter!(sp, gp1.g, gp1.p, label="Novel", markershape=:circle)
            scatter!(sp, gp0.g, gp0.p, label="Anceltral", markershape=:diamond)
            oname= @sprintf("%s_%.3d.pdf", basename, igen)
            Plots.pdf(oname)

            push!(gpdata,
                  (s.basename, epoch, igen,
                   mean(gp1.g), mean(gp1.p), std(gp1.g), std(gp1.p),
                   sum(var(geno1; dims=2)), sum(var(pheno1; dims=2)), mean(gp1.m1), mean(gp1.m0),
                   mean(gp0.g), mean(gp0.p), std(gp0.g), std(gp0.p),
                   sum(var(geno0; dims=2)), sum(var(pheno0; dims=2)), mean(gp0.m1), mean(gp0.m0))
                  )
        end
        sp = scatter(size=(600,600), 
                plot_title= s.basename,
                xlims=(-0.05, 1.05), ylims=(-0.05, 1.05),
                xticks=0:0.1:1, yticks=0:0.1:1,
                xlabel="Genotype", ylabel="Phenotype")
        scatter!(sp, gpdata.mg1, gpdata.mp1, xerror=gpdata.dg1, yerror=gpdata.dp1,
                 label="Novel", markershape=:circle, lc=:auto)
        scatter!(sp, gpdata.mg0, gpdata.mp0, xerror=gpdata.dg0, yerror=gpdata.dp0,
                 label="Anceltral", markershape=:diamond, lc=:auto)

        oname= @sprintf("%s_summary.pdf", basename)
        Plots.pdf(oname)
    end
    gpdata
end

function main()
    parsed_args = parse_commandline()
    ngen = parsed_args["ngen"]
    dir= parsed_args["outdir"]
    cname = parsed_args["out"]
    trajfiles = parsed_args["traj"]

    gpdata = mkdf()
    for traj in trajfiles
        append!(gpdata, proctraj(dir, ngen, traj))
    end

    CSV.write(cname, gpdata)
end

main()

