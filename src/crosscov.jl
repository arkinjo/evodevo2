using JLD2,CodecZlib
using ArgParse
using Statistics
using LinearAlgebra

include("EvoDevo2.jl")

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--gen"
        help = "generation number to be analyzed"
        arg_type = Int64
        default = 1
        
        "traj"
        help = "trajectory JLD2 file"
        arg_type = String
        required = true
        
    end

    return parse_args(s)
end

function anacc(label, cc, denv)
    nrow,ncol = size(cc)
    @printf("%s\tDims\t%d\t%d\n", label, nrow, ncol)
    F = svd(cc)
    
    stot = sum(F.S.^2) # total cross-covariance
    @printf("%s\tTOT\t%e\n", label, stot)
    ali = abs(F.U[:,1] ⋅ denv)
    @printf("%s\tAli\t%e\n", label, ali)

    for (i,s) in enumerate(F.S)
        @printf("%s\tSV\t%d\t%e\t%e\n", label, i, s, s^2/stot)
    end

    for (i,u) in enumerate(F.U[:,1])
        @printf("%s\tU1\t%d\t%e\n", label, i, u)
    end
end

function main()
    parsed_args = parse_commandline()
    igen = parsed_args["gen"]
    trajfile = parsed_args["traj"]
    jldopen(trajfile, "r") do file
        s = file["setting"]
        envs0 = file["envs0"]
        envs1 = file["envs1"]
        sel0 = selecting_envs(envs0,s)
        sel1 = selecting_envs(envs1,s)
        dsel = sel1 - sel0
        denvs = normalize(dsel)
                
        epoch = file["epoch"]
        pop0 = file[make_pop_name(Ancestral,igen)];
        pop1 = file[make_pop_name(Novel,igen)];

        pheno0 = get_selected_pheno_vecs(pop0, s)
        pheno1 = get_selected_pheno_vecs(pop1, s)

        geno0 = get_geno_vecs(pop0)
        geno1 = get_geno_vecs(pop1)

        cues0 = get_cue_vecs(pop0, s)
        cues1 = get_cue_vecs(pop1, s)

        pc0 = cov(pheno0, cues0; dims=2)
        pc1 = cov(pheno1, cues1; dims=2)

        anacc("PhenoCue_Nov", pc1, denvs)
        anacc("PhenoCue_Anc", pc0, denvs)

        pg0 = cov(pheno0, geno0; dims=2)
        pg1 = cov(pheno1, geno1; dims=2)

        anacc("PhenoGeno_Nov", pg1, denvs)
        anacc("PhenoGeno_Anc", pg0, denvs)
    end
end

main()

