using Printf
using LinearAlgebra
using ThreadsX

struct Population
    gen::Int64 # generation
    naenv::NAEnv
    indivs::Vector{Individual}
end

struct PopStats
    mismatch::Float64
    fitness::Float64
    ndev::Float64
    ppheno::Float64
    nparents::Int64
end

function PopStats(pop::Population, denv, env0, s::Setting)
    mismatch = 0.0
    fitness = 0.0
    ndev = 0.0
    ppheno = 0.0    
    npop = length(pop.indivs)
    npar = Dict()
    for indiv in pop.indivs
        mismatch += indiv.mismatch[]
        fitness += indiv.fitness[1]
        ndev += indiv.ndev[]
        f = get_selected_phenotype(indiv, s) - env0
        ppheno += dot(f, denv)
        npar[indiv.mom_id] = getkey(npar, indiv.mom_id, 0) + 1
    end
    mismatch /= npop
    fitness /= npop
    ndev /= npop
    ppheno /= npop
    PopStats(mismatch, fitness, ndev, ppheno, length(npar))
end

function Population(naenv::NAEnv, s::Setting)
    indivs = Vector{Individual}()
    for i = 1:s.max_pop
        push!(indivs, Individual(i, s))
    end
    Population(1, naenv, indivs)
end


function develop(pop::Population, envs::EnvironmentS, s::Setting)
    ThreadsX.foreach(pop.indivs) do indiv
        set_cues(indiv, envs, s)
        develop(indiv, envs, s)
    end
    set_relfit(pop)
end

function set_relfit(pop::Population)
    maxfit = mapreduce(ind -> ind.fitness[1], max, pop.indivs)
    for indiv in pop.indivs
        indiv.fitness[2] =  indiv.fitness[1]/maxfit
    end
end

function select(pop::Population, s::Setting)
    irange = 1:length(pop.indivs)
    parents = Vector{Individual}()
    npop = 0
    while npop <= s.max_pop
        indiv = rand(pop.indivs)
        if rand() < indiv.fitness[2]
            push!(parents, indiv)
            npop += 1
        end
    end
    parents
end

function reproduce(pop::Population, muts::Mutation, s::Setting) ::Vector{Individual}
    irange = 1:length(pop.indivs)
    parents = select(pop, s)
    dads = parents[1:2:end-1]
    moms = parents[2:2:end]
    ThreadsX.mapreduce(vcat, 1:length(dads), dads,moms) do i,dad,mom
        geno1, geno2 = mate(dad.genome, mom.genome)
        mutate(geno1, muts, s)
        mutate(geno2, muts, s)
        kid1 = Individual(2*i-1,dad.id,mom.id,geno1,s)
        kid2 = Individual(2*i,  mom.id,dad.id,geno2,s)
        [kid1,kid2]
    end
end


function evolve(mode, iepoch::Int64, ngen::Int64, pop1::Population,
                env0::EnvironmentS, env1::EnvironmentS,
                log, traj, s::Setting)
    indivs0 = deepcopy(pop1.indivs)
    pop0 = Population(1, Ancestral, indivs0)
    e0 = get_selecting_envs(env0, s)
    e1 = get_selecting_envs(env1, s)
    de = (e1 - e0)
    denv = de/dot(de, de)
    muts = Mutation(s)
    traj["envs0"] = env0
    traj["envs1"] = env1
    for igen = 1:ngen
        develop(pop1, env1, s)
        ps1 = PopStats(pop1, denv, e0, s)
        @printf(log, "%3d\t%3d", iepoch, igen)
        @printf(log, "\t%e\t%e\t%e\t%e\t%d", ps1.mismatch, ps1.fitness,
                ps1.ndev, ps1.ppheno, ps1.nparents)
        if mode == TestMode
            develop(pop0, env0, s)
            ps0 = PopStats(pop0, denv, e0, s)
            @printf(log, "\t%e\t%e\t%e\t%e\t%d",
                    ps0.mismatch, ps0.fitness, ps0.ndev, ps0.ppheno,
                    ps0.nparents)
            if traj != nothing
                name0 = @sprintf("pop0_%.3d", igen)
                name1 = @sprintf("pop1_%.3d", igen)
                traj[name0] = pop0
                traj[name1] = pop1
            end
            
            indivs0 = reproduce(pop1, muts, s)
            pop0 = Population(igen+1, Ancestral, indivs0)
            if igen % 10 == 0
                flush(log)
            end
        end
        @printf(log, "\n")
        indivs1 = reproduce(pop1, muts, s)
        pop1 = Population(igen+1, Novel, indivs1)
    end
    pop1
end

function train_epochs(nepoch::Int64, ngen::Int64, log, s::Setting)
    envs0 = make_environments(s)
    pop = Population(Novel, s)
    for iepoch = 1:nepoch
        envs1 = change_envS(envs0, iepoch + s.seed, s)
        pop = evolve(TrainMode, iepoch, ngen, pop, envs0, envs1,
                     log, nothing, s)
        envs0 = envs1
        flush(log)
    end
    envs0, pop
end

