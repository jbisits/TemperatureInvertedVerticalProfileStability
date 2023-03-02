using .VerticalProfileStability

sim_data_path = joinpath(SIM_DATADIR, "initial_ΔΘ_1.0_mu/")

files = glob("*.jld2", sim_data_path)
data_files = files[1:2]

## Get the time series of T, S, κ, buoyancy gradient and maximum density difference where
#  ΔΘ_thres is exceeded from the model output. By default ΔΘ_thres = 0.5, so for other
#  choices of initial ΔΘ this must be changed.

saved_timeseries = joinpath(sim_data_path, "output_timeseries.jld2")
save_timeseries!(saved_timeseries, data_files; ΔΘ_thres = 0.5)
