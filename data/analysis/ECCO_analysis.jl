using CairoMakie

ECCO_data = glob("*.nc", ECCO_datadir)
timestamps = Date(2007, 01, 01):Day(1):Date(2007, 01, 05)
series = RasterSeries(ECCO_data[1:5], Ti(timestamps); child = RasterStack)

ΔΘ_thres = 1.0
Δρ_max_series = series_max_Δρ(series, ΔΘ_thres)

Θₗ = series2vec(Δρ_max_series, :Θₗ)
get_lats(Δρ_max_series, :Θₗ) == get_lats(Δρ_max_series, :Δρ_static)
lats = get_lats(Δρ_max_series, :Θₗ)
Δρˢ = series2vec(Δρ_max_series, :Δρ_static)

fig = Figure()
ax = Axis(fig[1, 1];
          title = "Scatter of maximum density difference of a profile against\ntemperature of lower level where max dd was calculated",
          xlabel = "Θ (ᵒC) at lower level",
          xaxisposition = :top,
          ylabel = "Maximum Δρ (kgm⁻³) of profile")
#scatter!(ax, Θₗ, Δρˢ; markersize = 4) # no colour
sc = scatter!(ax, Θₗ, Δρˢ; markersize = 4, color = lats)
lines!(ax, zeros(4), 0:-1:-3; linestyle = :dash, linewidth = 3, color = :red)
lines!(ax, 0.5 .* ones(4), 0:-1:-3; linestyle = :dash, linewidth = 3, color = :orange)
lines!(ax, -1.88 .* ones(4), 0:-1:-3; linestyle = :dash, linewidth = 3, color = :green)
Colorbar(fig[1, 2], sc, label = "Latitude (ᵒN)")
fig
ylims!(ax, -0.5, maximum(Δρˢ)+0.01)
fig
