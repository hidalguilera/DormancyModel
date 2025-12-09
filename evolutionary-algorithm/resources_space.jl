using DrWatson
@quickactivate "dormancy"

using Agents, Agents.Schedulers, Distributions
import Agents: nextid, remove_agent!

mutable struct organism <: Agents.AbstractAgent

    id::Int
    pos::Dims{2}
    awake_rate::Float64
    awake_time::Float64
    death_time::Float64

end


mutable struct seme <: Agents.AbstractAgent

    id::Int
    pos::Dims{2}
    awake_rate::Float64
    awake_time::Float64
    death_time::Float64

end



Base.@kwdef mutable struct parameters_res_sim

    res::Int
    res_μ::Float64
    res_σ::Float64
    res_in::Vector{Float64}
    τ::Float64
    
    death_rate::Float64
    Δ::Float64
    p_mutate::Float64

    index::Int
    state::Int
    switch_time::Float64

    nseme::Int
    norganisms::Int

    time::Int

end



# Function to initialize the system       +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



function initialize_model(; res_μ = 100,  res_σ = 1, τ=10.0,
    n_organisms = 1000,  death_rate = 0.1, Δ = 0.1, p_mutate = 0.01,
    initial_a = 1
)

    space = Agents.GridSpace((1,1), periodic = true, metric=:euclidean)
    properties = parameters_res_sim(
        
        res_μ,
        res_μ,
        res_σ,
        zeros(2),
        τ,

        death_rate, Δ, p_mutate,

        1, -1, 0,

        n_organisms,
        100,

        0)

    model = Agents.ABM(Union{organism,seme}, space; properties, 
                      scheduler = Schedulers.Randomly(),
                      agent_step! = Organism_step!,
                      model_step! = model_step!)
    calcola_env_switch_time(0, model)
    define_resources(model)

    for id = 1:n_organisms
        death_time      = calcola_death_time(model.time, model)
        nagent = Agents.add_agent!(seme(id, (1,1), τ*initial_a, 0, death_time), model)
        calcola_awake_time(0, nagent)
        nagent.death_time      = calcola_death_time(nagent.awake_time, model)

    end
    
    return model

end



function define_resources(model::ABM)

    model.res       = Int(floor(rand(Uniform(model.res_μ-model.res_σ,model.res_μ+model.res_σ))))
    model.res_in    = [Int(floor(model.res_μ-model.res_σ)),Int(floor(model.res_μ+model.res_σ))]

end



function calcola_death_time(t, model)

    return t + log(1/rand())/(model.death_rate)

end



function calcola_awake_time(t, agent)

    agent.awake_time = t + log(1/rand())*agent.awake_rate#abs(rand(Normal(agent.awake_rate,1)))

end




# Functions for environmental dynamics       +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




function model_step!(model::ABM)
    
    model.time  += 1
    if model.time > model.switch_time
        model.state *= -1
        model.index = (model.state+1)/2+1
        model.switch_time = calcola_env_switch_time(model.time, model)
    end
    model.res   = model.res_in[model.index]
    model.nseme = count(a -> a isa seme, allagents(model))
    model.norganisms = count(a -> a isa organism, allagents(model))

end



function calcola_env_switch_time(t, model)

    model.switch_time =  t + log(1/rand())*model.τ

end


# Functions for individuals' dynamics       +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++





function Organism_step!(agent::organism, model::ABM; dt = 1)

    if agent.death_time < model.time
        print("Aiuto\n")
        if nagents(model)>1
            Agents.remove_agent!(agent, model)
        else
            agent.death_time += 2
        end

    elseif agent.death_time < (model.time+1)

            res_in = eat!(agent, model; dt=dt-(1+model.time - agent.death_time))
            n_offsprings = rand( Poisson( res_in) )

            for i in 1:n_offsprings
                reproduce!( agent, model)
            end
    
            if nagents(model)>1
                Agents.remove_agent!(agent, model)
            else
                agent.death_time += 2
            end
    
    else
        
        res_in = eat!(agent, model; dt=dt)
        n_offsprings = rand( Poisson( res_in) )

        for i in 1:n_offsprings
            reproduce!( agent, model)
        end
    
    end
end



function Organism_step!(agent::seme, model::ABM)

    if agent.awake_time < model.time + 1
        
        if agent.death_time > (model.time+1)# && conta_svegli(allagents(model)) < model.res_μ
            new_id = nextid(model)
            Agents.add_agent!(organism(new_id, (1,1), 
                agent.awake_rate, 
                agent.awake_time, 
                agent.death_time), 
                model)

            Organism_step!(model[new_id], model; dt=1+model.time - agent.awake_time)

        end

        if nagents(model)>1
            Agents.remove_agent!(agent, model)
        else
            agent.death_time += 2
        end

    end

end


function eat!(agent::organism, model::ABM; dt = 1.0)

    res_in            = rand( Poisson( dt * model.res / (model.norganisms+1)) )
    return res_in

end



function reproduce!(agent::organism, model::ABM)

    new_awake_rate          = mutate_awake_rate(agent.awake_rate, model)   

    nagent = Agents.add_agent!(seme(nextid(model), (1,1), 
                new_awake_rate, 
                model.time, 
                model.time), 
                model)

    calcola_awake_time(model.time + 1 - rand(), nagent)
    nagent.death_time = calcola_death_time(nagent.awake_time, model)
    
    if nagent.awake_time < model.time + 1
        Organism_step!(nagent, model)
    end
    
    return

end


function mutate_awake_rate(awake_rate::Float64, model::ABM)
    
    new_awake_rate     = copy(awake_rate)
    
    if rand() < model.p_mutate
        new_awake_rate    *= rand(Normal(1.00,model.Δ))
        new_awake_rate     = (new_awake_rate>0 ? new_awake_rate : awake_rate)
    end

    return new_awake_rate

end


