using .VerticalProfileStability
using JLD2, Statistics

const ARGO_OUTPUT = jldopen(joinpath(@__DIR__, "ARGO_extracted.jld2"))

## Temperature inverted plots
# set ylimits, much faster to extract then plot data then plot and use lims! on axis.
xlimits = (-1.88, 10)
ylimits = (-0.1, 0.1)
# Density thresholds
Θₗ_range = range(-1.8, 5; length = 100)
ΔΘ_thres = [0.5, 1.0, 2.0, 3.0]
## Same axis, thresholds in different colours
fig = Figure(size = (700, 700))
ax = Axis(fig[1, 1];
        xlabel = "Θ (ᵒC)",
        xaxisposition = :top,
        title = "Maximum static density difference between two vertically\nspaced levels of a profile against temperature of lower level",
        ylabel = "Δρ (kgm⁻³)")
ylims!(ax, ylimits)
colors = [:blue, :orange, :red, :green]
for (i, key) ∈ enumerate(keys(ARGO_OUTPUT))

    Θᵤ = collect(skipmissing(ARGO_OUTPUT[key]["Θᵤ"]))
    Θₗ = collect(skipmissing(ARGO_OUTPUT[key]["Θₗ"]))
    find_inverted = findall(Θᵤ .< Θₗ)
    Δρˢ = collect(skipmissing(ARGO_OUTPUT[key]["Δρˢ"]))
    sc = scatter!(ax, Θₗ[find_inverted], Δρˢ[find_inverted];
                  color = colors[i], markersize = 5)

    Sₗ_mean = mean(collect(skipmissing(ARGO_OUTPUT[key]["Sₗ"][find_inverted])))
    pₘ = 0.5 .* (collect(skipmissing(ARGO_OUTPUT[key]["pₗ"][find_inverted])) .+
                 collect(skipmissing(ARGO_OUTPUT[key]["pᵤ"][find_inverted])))
    pₘ_mean = mean(pₘ)
    αₗ = gsw_alpha.(Sₗ_mean, Θₗ_range, pₘ_mean)
    βₗ = gsw_beta.(Sₗ_mean, Θₗ_range, pₘ_mean)
    Δρ_thres = @. gsw_rho.(Sₗ_mean - (αₗ / βₗ) * ΔΘ_thres[i], Θₗ_range - ΔΘ_thres[i],
                          pₘ_mean) - gsw_rho.(Sₗ_mean, Θₗ_range, pₘ_mean)
    lines!(ax, Θₗ_range, Δρ_thres; color = colors[i],
           label = "ΔΘ = $(ΔΘ_thres[i])ᵒC")

end
#axislegend(ax; position = :rb)
Legend(fig[2, 1], ax, "Δρ threshold for", orientation = :horizontal)
fig
save(joinpath(PLOTDIR, "ARGO", "ΔΘ_thres_all.png"), fig)
ylims!(ax, -0.1, 0.01)
fig
save(joinpath(PLOTDIR, "ARGO", "ΔΘ_thres_all_zoom.png"), fig)
## Full plot
# Temperature colourbar is wrong here, need a rethink for this and the pressure difference
keys_mat = reshape(keys(ARGO_OUTPUT), 2, 2)
colourbar_vars = ("lats", "Δp", "ΔΘ")
colourbar_col = (:viridis, :batlow, :thermal)
# These have been found by looking at the subset for each key in the dictionary and taking
# the `extrema`. The code is not here though.
colourbar_ranges = ((-80, -60), (20, 410), (-3.5, -0.5))
choose_var = 2
## Setup axis
fig = Figure(size = (1400, 1400))
ax = [Axis(fig[i, j];
          xlabel = "Θ (ᵒC)",
          xaxisposition = :top,
          ylabel = "Δρ (kgm⁻³)",
          title = keys_mat[i, j]*"ᵒC")
      for i ∈ 1:2, j ∈ 1:2]
for i ∈ 1:2
    linkxaxes!(ax[i, :]...)
    linkyaxes!(ax[:, i]...)
    hidexdecorations!(ax[2, i], grid = false)
    hideydecorations!(ax[i, 2], grid = false)
end
for i ∈ 1:2, j ∈ 1:2
    ylims!(ax[i, j], ylimits)
end
fig

## Plot
for (i, key) ∈ enumerate(keys_mat)

    Θᵤ = collect(skipmissing(ARGO_OUTPUT[key]["Θᵤ"]))
    Θₗ = collect(skipmissing(ARGO_OUTPUT[key]["Θₗ"]))
    Δρˢ = collect(skipmissing(ARGO_OUTPUT[key]["Δρˢ"]))
    plot_var = round.(ARGO_OUTPUT[key][colourbar_vars[choose_var]][find]; digits = 1)
    find_inverted = findall(Θᵤ .< Θₗ)
    # plot_var = round.(ARGO_OUTPUT[key][colourbar_vars[choose_var]][find_inverted]; digits = 1)
    cbar_lab =  if colourbar_vars[choose_var] == "lats"
                   "Latitude (°N)"
                elseif colourbar_vars[choose_var] == "Δp"
                    "Δp (dbar)"
                elseif colourbar_vars[choose_var] == "ΔΘ"
                    "ΔΘ (°C)"
                end

    sc = scatter!(ax[i], Θₗ[find_inverted], Δρˢ[find_inverted];
                  markersize = 4,
                  #color = plot_var,
                  color = lats[find_inverted],
                  colormap = colourbar_col[choose_var],
                  colorrange = colourbar_ranges[choose_var])
    # lines!(ax[i], Θ_lower_range, Δρ_thres;
    #        color = colors[i],
    #        label = "ΔΘ = $(ΔΘ_range[1])ᵒC")

    if i == 4
        Colorbar(fig[:, 3], sc, label = cbar_lab)
    end
end
fig
# lin_elements = [LineElement(color = col) for col ∈ colors]
# leg_labels = "ΔΘ = " .* string.([0.5, 1.0, 2.0, 3.0]) .* "°C"
# Legend(fig[3, :], lin_elements, leg_labels, "Δρ threshold for", orientation = :horizontal)
close(ARGO_OUTPUT)
