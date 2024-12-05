using ArgParse

include("EvoDevo2.jl")

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--seed"
        help = "random seed (13579)"
        arg_type = Int64
        default = 13579

        "traj"
        nargs = '+'
        help = "trajectory JDL2 file"
        arg_type = String
        required = true
    end

    return parse_args(s)
end

function get_cov(pop0, pop1, envs0, envs1, s)
    initialize(pop0, s)
    initialize(pop1, s)

    develop(pop1, envs0, s)
    develop(pop1, envs1, s)

    pheno0 = get_selected_pheno_vecs(pop0, s)
    pheno1 = get_selected_pheno_vecs(pop1, s)
    dpheno = map((p0,p1) -> p1 - p0, pheno0, pheno1)
    cov0 = norm(cov(pheno0, pheno0; dims=2))
    cov1 = norm(cov(pheno1, pheno1; dims=2))
    dcov = norm(cov(dpheno, dpheno; dims=2))
    return cov0,cov1,dcov
end

function main()
    parsed_args = parse_commandline()
    seed = parsed_args["seed"]
    trajfiles = parsed_args["traj"]
    for file in trajfiles
        s, envs0, envs1, pop0, pop1 =
            jldopen(resfile, "r") do file
                file["setting"],
                file["envs0"],
                file["envs1"],
                file[make_pop_name(Ancestral,1)],
                file[make_pop_name(Novel,1)]
            end

        env_noise = s.env_noise
        s.env_noise = 0
        varG0,varG1,varD0 = get_cov(pop0, pop1, envs0, envs1, s)
        s.env_noise = 0.05
        varP0,varP1,varD1 = get_cov(pop0, pop1, envs0, envs1, s)

        println("heritability=", varD0/varD1, varG0/varP0, varG1/varP1)
    end
end

main()

