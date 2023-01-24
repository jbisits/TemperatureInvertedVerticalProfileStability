using CairoMakie, JLD2
## Choose which simulation
## Initial temperature -1.85 and 0.5ᵒC
sim_output = joinpath(sim_datadir, "initial_Θ_minus1_85_0_5")
## Initial temperature 18.15 and 20.5ᵒC
sim_output = joinpath(sim_datadir, "initial_Θ_18_15_20_5")

## Open chosen simulation
saved_ts = jldopen(joinpath(sim_output, "output_timeseries.jld2"))
t, S_ts, T_ts, κ_ts, Δρ_ts = saved_ts["t"], saved_ts["S_ts"], saved_ts["T_ts"],
                             saved_ts["κ_ts"], saved_ts["Δρ_ts"]
close(saved_ts)
z = -500:5:-5
length(z)
## Initial conditions
S₀, T₀ = [S_ts[:, 1, i] for i ∈ 1:length(S_ts[1, 1, :])],
         [T_ts[:, 1, i] for i ∈ 1:length(T_ts[1, 1, :])]
lines(T₀[1], z)
lines(S₀[1], z)
## Diffusivity time series
heatmap(t, z, κ_ts[:, :, 1]')
## Δρ static and cabbeling time series
lines(t, Δρ_ts["Δρ_c"][:, 1])
