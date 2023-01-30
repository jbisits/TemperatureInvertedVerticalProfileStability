## Simulation with no diffusion for testing, no need to include in the test as data is
#  saved.
# savepath = joinpath(@__DIR__, "test_sim_data")
# convective_κ, background_κ = 0.0, 0.0
# num_ics = 2
# run_model(; convective_κ, background_κ, savepath, num_ics)

# Test
sim_data_path = joinpath(@__DIR__, "test_sim_data")
saved_sim_data = glob("*.jld2", sim_data_path)[1:2]
T_test = Array{Vector}(undef, 2, 2)
S_test = Array{Vector}(undef, 2, 2)
for (i, data) ∈ enumerate(saved_sim_data)

    temp_output = jldopen(data)
    # Initial and final temperature
    T_test[i, 1] = reshape(temp_output["timeseries"]["T"][t[1]], :)
    T_test[i, 2] = reshape(temp_output["timeseries"]["T"][t[end]], :)
    # Initial and final salt
    S_test[i, 1] = reshape(temp_output["timeseries"]["S"][t[1]], :)
    S_test[i, 2] = reshape(temp_output["timeseries"]["S"][t[end]], :)
    close(temp_output)

end
