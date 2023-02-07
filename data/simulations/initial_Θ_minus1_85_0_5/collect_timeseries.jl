using JLD2
using .VerticalProfileStability

sim_data_path = joinpath(sim_datadir, "initial_Θ_minus1_85_0_5/")

files = glob("*.jld2", sim_data_path)
data_files = files[1:11]
data_files = order_files(data_files)

## Get the time series of T, S, κ, buoyancy gradient and maximum density difference where
#  ΔΘ_thres is exceeded from the model output. By default ΔΘ_thres = 0.5 but will look at
#  other values.

saved_timeseries = joinpath(sim_data_path, "output_timeseries.jld2")
save_timeseries!(saved_timeseries, data_files)

## ΔΘ_thres = 0.25
ΔΘ_thres_0_25_timeseries = joinpath(sim_data_path, "ΔΘ_thres_0_25_timeseries.jld2")
save_Δρ_timeseries!(ΔΘ_thres_0_25_timeseries, data_files; ΔΘ_thres = 0.25)

## ΔΘ_thres = 1.0
ΔΘ_thres_1_timeseries = joinpath(sim_data_path, "ΔΘ_thres_1_timeseries.jld2")
save_Δρ_timeseries!(ΔΘ_thres_1_timeseries, data_files; ΔΘ_thres = 1.0)

## ΔΘ_thres = 2.0
ΔΘ_thres_2_timeseries = joinpath(sim_data_path, "ΔΘ_thres_2_timeseries.jld2")
save_Δρ_timeseries!(ΔΘ_thres_2_timeseries, data_files; ΔΘ_thres = 2.0)
