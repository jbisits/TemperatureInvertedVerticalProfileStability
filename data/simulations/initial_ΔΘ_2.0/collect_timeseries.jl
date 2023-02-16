using JLD2
using .VerticalProfileStability

sim_data_path = joinpath(sim_datadir, "initial_ΔΘ_2.0/")

files = glob("*.jld2", sim_data_path)
data_files = files[1:3]

## Get the time series of T, S, κ, buoyancy gradient and maximum density difference where
#  ΔΘ_thres is exceeded from the model output. By default ΔΘ_thres = 0.5, so for other
#  choices of initial ΔΘ this must be changed.

saved_timeseries = joinpath(sim_data_path, "output_timeseries.jld2")
save_timeseries!(saved_timeseries, data_files; ΔΘ_thres = 1.0)
