using ArgParse
include("EvoDevo2.jl")

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--out"
        default = "cc.csv"
        help = "output CSV filename"
        arg_type = String
        
        "--traj"
        nargs = '+'
        help = "trajectory JLD2 files"
        arg_type = String
    end

    return parse_args(s)
end

function anacc(label, epoch, novanc, geno, cue, denv, pprj)
    G = svd(geno)
    gtot = sum(G.S.^2)
    gali = abs(G.U[:,1] ⋅ denv)
    gs1 = G.S[1]

    C = svd(cue)
    ctot = sum(C.S.^2)
    cali = abs(C.U[:,1] ⋅ denv)
    cs1 = C.S[1]

    gc = G.U[:,1] ⋅ C.U[:,1]

    pave,pstd = pprj
    (label, epoch, novanc,
     pave, pstd,
     gtot, gs1, (gs1^2)/gtot, gali,
     ctot, cs1, (cs1^2)/ctot, cali,
     gc)
end

function proctraj(trajfile)
    println(stderr, "# $trajfile")

    s,epoch,envs0,envs1,pop0,pop1 =
        jldopen(trajfile, "r") do file
            file["setting"],
            file["epoch"],
            file["envs0"],
            file["envs1"],
            file[make_pop_name(Ancestral,1)],
            file[make_pop_name(Novel,1)]
        end

    sel0 = selecting_envs(envs0,s)
    sel1 = selecting_envs(envs1,s)
    dsel = sel1 - sel0
    paxis = dsel/(dsel ⋅ dsel)
    denvs = normalize(dsel)

    pheno0 = get_selected_pheno_vecs(pop0, s)
    pheno1 = get_selected_pheno_vecs(pop1, s)

    function project(phenos)
        n = size(phenos)[2]
        prj = ThreadsX.map(i -> (phenos[:,i] - sel0) ⋅ paxis, 1:n)
        Statistics.mean(prj),Statistics.std(prj)
    end
    pprj0 = project(pheno0)
    pprj1 = project(pheno1)

    geno0 = get_geno_vecs(pop0)
    geno1 = get_geno_vecs(pop1)

    cues0 = get_cue_vecs(pop0, s)
    cues1 = get_cue_vecs(pop1, s)

    pc0 = cov(pheno0, cues0; dims=2)
    pg0 = cov(pheno0, geno0; dims=2)
    anc = anacc(s.basename, epoch, "Anc", pg0, pc0, denvs, pprj0)
    
    pc1 = cov(pheno1, cues1; dims=2)
    pg1 = cov(pheno1, geno1; dims=2)
    nov = anacc(s.basename, epoch, "Nov", pg1, pc1, denvs, pprj1)

    (anc, nov)
end


function main()
    parsed_args = parse_commandline()
    outfile = parsed_args["out"]
    trajfiles = parsed_args["traj"]
    data = DataFrame(basename = String[],
                     epoch=Int64[],
                     novanc = String[], # Nov or Anc

                     pave=Float64[], # mean projected phenotype
                     pstd=Float64[], # stdev of projected phenotype
                     
                     tot_geno=Float64[], # total cross-cov
                     s1_geno=Float64[], # first singular value
                     ps1_geno=Float64[], # % first singular value
                     ali_geno=Float64[], # alignment vs. denv

                     tot_cue=Float64[], # total cross-cov
                     s1_cue=Float64[], # first singular value
                     ps1_cue=Float64[], # % first singular value
                     ali_cue=Float64[], # alignment vs. denv

                     geno_cue=Float64[] # U1(geno) \cdot U1(cue)
                     )

    for (anc,nov) in ThreadsX.map(proctraj, trajfiles)
        push!(data, anc)
        push!(data, nov)
    end
    CSV.write(outfile, data)

end

main()

