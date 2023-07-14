using JLD2,CodecZlib
include("EvoDevo2.jl")

s = default_setting("full", 5, 13)
s.max_pop = 1000
s.num_cell_x = 2
s.num_cell_y = 2
set_omegas(s)

envs0 = make_environments(s);
envs1 = change_envS(envs0, 111 + s.seed, s);

if ARGS[1] == "no" 
    println("### No threads ###")
    pop = Population(Novel, s);

    for i = 1:10
        @time for indiv in pop.indivs
            set_cues(indiv, envs1, s)
            develop(indiv, envs1, s)
        end
    end
else
    println("### With threads ###")
    pop = Population(Novel, s);
    for i = 1:10
        @time Threads.@threads for indiv in pop.indivs
            set_cues(indiv, envs1, s)
            develop(indiv, envs1, s)
        end
    end
end

a = ones(200);
b = map(a[1:100], a[101:200]) do x,y
    x + y
end
