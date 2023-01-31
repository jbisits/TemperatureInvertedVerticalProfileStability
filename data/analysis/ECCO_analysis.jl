using CairoMakie
using .VerticalProfileStability

ECCO_data = glob("*.nc", ECCO_datadir)
timestamps = Date(2007, 01, 01):Day(1):Date(2007, 12, 31)
series = RasterSeries(ECCO_data, Ti(timestamps); child = RasterStack)

## Compute max denisty difference for given ΔΘ_thres and extract data to plot
ΔΘ_thres = 2.0
Δρ_max_series = series_max_Δρ(series, ΔΘ_thres)

Θₗ = series2vec(Δρ_max_series, :Θₗ)
lats = get_lats(Δρ_max_series, :Θₗ)
Δρˢ = series2vec(Δρ_max_series, :Δρ_static)
Δρᶜ = series2vec(Δρ_max_series, :Δρ_cab)
Δρ = [Δρˢ, Δρᶜ]

## Full plot
fig = Figure(size = (1200, 600))
dd = ["static", "cabbeling"]
ax = [Axis(fig[1, i];
          title = "Maximum $(dd[i]) density difference of a\nprofile against temperature of lower level\nwhere max dd was calculated with ΔΘ = $ΔΘ_thres",
          subtitle = "ECCO data January 2007",
          xlabel = "Θ (ᵒC) at lower level",
          xaxisposition = :top,
          ylabel = "Maximum Δρ (kgm⁻³) of profile") for i ∈ 1:2]
linkyaxes!(ax...)
hideydecorations!(ax[2])
for (i, Δρ_) ∈ enumerate(Δρ)
    sc = scatter!(ax[i], Θₗ, Δρ_; markersize = 4, color = lats)
    #scatter!(ax[i], Θₗ, Δρ; markersize = 4) # no colour
    lines!(ax[i], zeros(4), 0:-1:-3;
           linestyle = :dash, linewidth = 3, color = :red)
    lines!(ax[i], 0.5 .* ones(4), 0:-1:-3;
           linestyle = :dash, linewidth = 3, color = :orange)
    lines!(ax[i], -1.88 .* ones(4), 0:-1:-3;
           linestyle = :dash, linewidth = 3, color = :green)
end
Colorbar(fig[2, :], sc, label = "Latitude (ᵒN)", vertical = false, flipaxis = false)
#fig
save(joinpath(plotdir, "ECCO", "Jan2007_ΔΘ_thres_2.png"), fig)
ylims!(ax[1], -0.5, 0)
ylims!(ax[2], -0.5, 0)
#fig
save(joinpath(plotdir, "ECCO", "Jan2007_ΔΘ_thres_2_zoom.png"), fig)
ylims!(ax[1], -0.1, 0)
ylims!(ax[2], -0.1, 0)
#fig
save(joinpath(plotdir, "ECCO", "Jan2007_ΔΘ_thres_2_zoomzoom.png"), fig)
