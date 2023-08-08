using JLD2,CodecZlib
include("EvoDevo2.jl")
using ThreadsX

s = default_setting("full", 5, 13)
s.max_pop = 1000
s.num_cell_x = 2
s.num_cell_y = 2
set_omegas(s)

envs0 = make_environments(s);
envs1 = change_envS(envs0, s);

if ARGS[1] == "no" 
    println("### No threads ###")
    pop = Population(Novel, s);

    for i = 1:10
        @time for indiv in pop.indivs
            set_cues(indiv, envs1, s)
            develop(indiv, envs1, s)
        end
    end
elseif ARGS[1] == "t"
    println("### With @threads ###")
    pop = Population(Novel, s);
    for i = 1:10
        @time Threads.@threads for indiv in pop.indivs
            set_cues(indiv, envs1, s)
            develop(indiv, envs1, s)
        end
    end
else
    println("### With ThreadsX ###")
    pop = Population(Novel, s);
    for i = 1:10
        @time ThreadsX.foreach(pop.indivs) do indiv
            set_cues(indiv, envs1, s)
            develop(indiv, envs1, s)
        end
    end
end

