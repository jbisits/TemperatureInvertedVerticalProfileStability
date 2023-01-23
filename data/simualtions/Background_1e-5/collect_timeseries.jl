using .VerticalProfileStability

sim_data_path = joinpath(sim_datadir, "Background_1e-5_withdiff")

files = glob("*.jld2", sim_data_path)
data_files = files[1:11]
data_files = order_files(data_files)

## Get the density difference where ΔΘ_thres is exceeded.

saved_densitydiff = joinpath(sim_data_path, "density_diff_thetathres.jld2")
save_densitydiff!(saved_densitydiff, data_files)

## Get the time series of T, S, κ and buoyancy from the model output

saved_timeseries = joinpath(sim_data_path, "output_timeseries.jld2")
save_timeseries!(saved_timeseries, data_files)
