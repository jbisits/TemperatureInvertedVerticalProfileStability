using JLD2
using .VerticalProfileStability

sim_data_path = joinpath(sim_datadir, "initial_Θ_minus1_85_0_5/")

files = glob("*.jld2", sim_data_path)
data_files = files[1:11]
data_files = order_files(data_files)

## Get the time series of T, S, κ, buoyancy gradient and maximum density difference where
#  ΔΘ_thres is exceeded from the model output.

saved_timeseries = joinpath(sim_data_path, "output_timeseries.jld2")
save_timeseries!(saved_timeseries, data_files)
