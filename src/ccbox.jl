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

function anacc(epoch, novanc, geno, cue, denv)
    G = svd(geno)
    gtot = sum(G.S.^2)
    gali = abs(G.U[:,1] ⋅ denv)
    gs1 = G.S[1]

    C = svd(cue)
    ctot = sum(C.S.^2)
    cali = abs(C.U[:,1] ⋅ denv)
    cs1 = C.S[1]

    gc = G.U[:,1] ⋅ C.U[:,1]
    (epoch, novanc,
     gtot, gs1, (gs1^2)/gtot, gali,
     ctot, cs1, (cs1^2)/ctot, cali,
     gc)
end

function proctraj(trajfile)
    jldopen(trajfile, "r") do file
        s = file["setting"]
        envs0 = file["envs0"]
        envs1 = file["envs1"]
        sel0 = selecting_envs(envs0,s)
        sel1 = selecting_envs(envs1,s)
        denvs = normalize(sel1 - sel0)

        epoch = file["epoch"]
        pop0 = file[make_pop_name(Ancestral,1)];
        pop1 = file[make_pop_name(Novel,1)];

        pheno0 = get_selected_pheno_vecs(pop0, s)
        pheno1 = get_selected_pheno_vecs(pop1, s)

        geno0 = get_geno_vecs(pop0)
        geno1 = get_geno_vecs(pop1)

        cues0 = get_cue_vecs(pop0, s)
        cues1 = get_cue_vecs(pop1, s)

        pc0 = cov(pheno0, cues0; dims=2)
        pg0 = cov(pheno0, geno0; dims=2)
        anc = anacc(epoch, "Anc", pg0, pc0, denvs)
        
        pc1 = cov(pheno1, cues1; dims=2)
        pg1 = cov(pheno1, geno1; dims=2)
        nov = anacc(epoch, "Nov", pg1, pc1, denvs)

        (anc, nov)
    end

end


function main()
    parsed_args = parse_commandline()
    outfile = parsed_args["out"]
    trajfiles = parsed_args["traj"]
    data = DataFrame(epoch=Int64[],
                     novanc = String[], # Nov or Anc
                     
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

    dlst =
        ThreadsX.map(trajfiles) do trajfile
            proctraj(trajfile)
        end
    for (anc,nov) in dlst
        push!(data, anc)
        push!(data, nov)
    end
    CSV.write(outfile, data)

end

main()

