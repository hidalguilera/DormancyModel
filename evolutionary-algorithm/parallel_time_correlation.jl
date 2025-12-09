using Distributed

# Add workers before loading packages
if nprocs() == 1
    addprocs(16; exeflags="--project=.")  # Add 16 worker processes with project environment
end

# Load packages on all processes
@everywhere begin
    using Pkg
    Pkg.activate(".")  # Activate current project
    
    using Plots, ProgressBars, DataFrames, StatsPlots, Dates, CSV
    using Agents, Agents.Schedulers, Distributions
    import Agents: nextid, remove_agent!
    
    # Include the simulation code
    include(joinpath(dirname(@__DIR__), "src", "resources_space.jl"))
    
    # Define plotting directory function
    plotsdir(args...) = joinpath(dirname(@__DIR__), "plots", args...)
    datadir(args...) = joinpath(dirname(@__DIR__), "data", args...)
end

# Define functions on all processes
@everywhere begin
    steps = 200000
    
    avg_awake_rate_organisms(model) = begin
        organisms = [agent for agent in allagents(model) if agent isa organism]
        isempty(organisms) ? 0.0 : mean([agent.awake_rate for agent in organisms])
    end
    
    conta_svegli(as) = sum([a isa organism for a in as])
    awake_rate(a) = a.awake_rate
    
    mdata = [(nagents), (avg_awake_rate_organisms)]
    adata = [(awake_rate, mean)]
    
    # Function to run a single simulation with replicate ID
    function run_single_simulation(sigma, initial_a, replicate_id; steps=200000)
        try
            println("Worker $(myid()): Running simulation with sigma=$sigma, initial_a=$initial_a, replicate=$replicate_id")
            
            m = initialize_model(; res_μ = 100, res_σ = sigma, τ=4.0,
                n_organisms = 1000, death_rate = 0.1, Δ = 0.1, p_mutate = 0.01,
                initial_a = initial_a
            )

            dataa, datam = run!(m, steps; adata=adata, mdata=mdata)
            
            # Calculate final awake rate (mean of last 1000 steps)
            final_awake_rate = mean(dataa[!,:mean_awake_rate][max(1, end-1000):end])
            
            return (sigma = sigma, 
                   init_a = dataa[!,:mean_awake_rate][1], 
                   final_awake_rate = final_awake_rate,
                   replicate = replicate_id,
                   worker_id = myid())
        catch e
            println("Worker $(myid()): Error in simulation with sigma=$sigma, initial_a=$initial_a, replicate=$replicate_id: $e")
            return (sigma = sigma, 
                   init_a = initial_a, 
                   final_awake_rate = NaN,
                   replicate = replicate_id,
                   worker_id = myid())
        end
    end
end

# Main execution on master process
println("Starting parallel simulations on $(nprocs()) processes ($(nworkers()) workers)")

# Create parameter combinations
sigma_values = [90.0]
initial_a_values = 10.0 .^ LinRange(-1.95, 0.75, 17)
n_replicates = 10

# Create all parameter combinations with replicates
param_combinations = [(sigma, initial_a, replicate) 
                     for sigma in sigma_values 
                     for initial_a in initial_a_values 
                     for replicate in 1:n_replicates]

println("Total number of simulations: $(length(param_combinations))")
println("Parameter combinations: $(length(sigma_values)) σ × $(length(initial_a_values)) initial_a × $n_replicates replicates")
println("Simulations per worker: ~$(ceil(length(param_combinations) / nworkers()))")

# Run simulations in parallel
println("Starting parallel execution...")
start_time = time()

# Use pmap for automatic load balancing
results_tuples = pmap(param_combinations) do (sigma, initial_a, replicate)
    run_single_simulation(sigma, initial_a, replicate; steps=steps)
end

end_time = time()
elapsed_time = end_time - start_time

println("All simulations completed!")
println("Total execution time: $(round(elapsed_time, digits=2)) seconds")
println("Average time per simulation: $(round(elapsed_time / length(param_combinations), digits=2)) seconds")

# Convert results to DataFrame
results = DataFrame(results_tuples)

# Remove failed simulations (NaN results)
valid_results = filter(row -> !isnan(row.final_awake_rate), results)
failed_count = nrow(results) - nrow(valid_results)

if failed_count > 0
    println("Warning: $failed_count simulations failed and were excluded")
end

println("Successfully completed simulations: $(nrow(valid_results))")

# Analysis and plotting (updated for replicates)
println("Creating plots...")

# Group results for analysis - now we have replicates to aggregate
# First, let's get statistics for each parameter combination across replicates
replicate_stats = combine(groupby(valid_results, [:sigma, :init_a]), 
    :final_awake_rate => mean => :final_awake_rate_mean,
    :final_awake_rate => std => :final_awake_rate_std,
    :final_awake_rate => (x -> length(x)) => :n_replicates)

# Also save individual replicate results for detailed analysis
println("Replicate statistics:")
println("  Mean replicates per parameter combination: $(mean(replicate_stats.n_replicates))")
println("  Min replicates: $(minimum(replicate_stats.n_replicates))")
println("  Max replicates: $(maximum(replicate_stats.n_replicates))")

# Save results - both individual replicates and aggregated statistics
timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")

# Save individual replicate results
results_file_individual = datadir("parallel_time_correlation_individual_replicates_$(timestamp).csv")
mkpath(dirname(results_file_individual))
CSV.write(results_file_individual, valid_results)
println("Individual replicate results saved to: $results_file_individual")

# Save aggregated statistics
results_file_stats = datadir("parallel_time_correlation_replicate_stats_$(timestamp).csv") 
CSV.write(results_file_stats, replicate_stats)
println("Replicate statistics saved to: $results_file_stats")

# Summary statistics
println("\nSummary Statistics:")
println("Parameter ranges:")
println("  Sigma: $(minimum(valid_results.sigma)) to $(maximum(valid_results.sigma))")
println("  Initial awake rate: $(minimum(valid_results.init_a)) to $(maximum(valid_results.init_a))")
println("  Final awake rate: $(minimum(valid_results.final_awake_rate)) to $(maximum(valid_results.final_awake_rate))")
println("Replication summary:")
println("  Total individual simulations: $(nrow(valid_results))")
println("  Number of parameter combinations: $(nrow(replicate_stats))")
println("  Replicates per combination: $n_replicates (target)")
println("  Mean final awake rate across all: $(round(mean(valid_results.final_awake_rate), digits=3))")
println("  Std final awake rate across all: $(round(std(valid_results.final_awake_rate), digits=3))")

# Worker usage summary
worker_usage = combine(groupby(valid_results, :worker_id), nrow => :simulations_completed)
println("\nWork distribution:")
for row in eachrow(worker_usage)
    println("  Worker $(row.worker_id): $(row.simulations_completed) simulations")
end

# Cleanup
println("\nRemoving worker processes...")
rmprocs(workers())
println("Done!")