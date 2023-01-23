using CairoMakie
using .VerticalProfileStability
using .VerticalProfileStability.MaximumDensityDifference
import .VerticalProfileStability.MaximumDensityDifference: series_max_Δρ

ECCO_data = glob("*.nc", ECCO_datadir)
timestamps = Date(2007, 01, 01):Day(1):Date(2007, 01, 02)
series = RasterSeries(ECCO_data[1:2], Ti(timestamps); child = RasterStack)
dims(series[1])

ΔΘ_thres = 1.0
Δρ_max_series = series_max_Δρ(series, ΔΘ_thres)

Θₗ = collect(skipmissing(Δρ_max_series[1][:Θₗ]))
Δρˢ = collect(skipmissing(Δρ_max_series[1][:Δρ_static]))

fig = Figure()
ax = Axis(fig[1, 1];
          title = "Scatter of maximum density difference of a profile against\ntemperature of lower level where max dd was calculated",
          xlabel = "Θ (ᵒC) at lower level",
          xaxisposition = :top,
          ylabel = "Maximum Δρ (kgm⁻³) of profile")
scatter!(ax, Θₗ, Δρˢ; markersize = 4)
lines!(ax, zeros(4), 0:-1:-3; linestyle = :dash, linewidth = 3, color = :red)
lines!(ax, 1 .* ones(4), 0:-1:-3; linestyle = :dash, linewidth = 3, color = :orange)
fig
ylims!(ax, -0.1, maximum(Δρˢ))
fig
