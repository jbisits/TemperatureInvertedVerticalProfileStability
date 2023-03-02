using .VerticalProfileStability
using Statistics

ARGO_OUTPUT = jldopen(joinpath(@__DIR__, "ARGO_extracted.jld2"))

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
colours = [:blue, :orange, :red, :green]
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

## Histogram
# This is lazy as I have just used the mean Δρ_thres, better would be to bin by temperature
# and count how many exceed the density threshold in a given bin? Have done something
# similar in the `ECCO_invertedΔΘ_extracted.jl` script.
using StatsBase, LinearAlgebra

ΔΘ_keys = keys(ARGO_OUTPUT)
bin_width = 0.01

fig = Figure(size = (1200, 1200))
ax = [Axis(fig[i, j];
        xlabel = "Δρ (kgm⁻³)"
        )
      for i ∈ 1:2, j ∈ 1:2]
less_thres = Vector{Float64}(undef, 4)
over_thres = Vector{Float64}(undef, 4)

for (i, key) ∈ enumerate(ΔΘ_keys)
    Θₗ = collect(skipmissing(ARGO_OUTPUT[key]["Θₗ"]))
    Θᵤ = collect(skipmissing(ARGO_OUTPUT[key]["Θᵤ"]))
    find_inverted = findall(Θᵤ .< Θₗ)
    Δρˢ = collect(skipmissing(ARGO_OUTPUT[key]["Δρˢ"][find_inverted]))
    Sₗ_mean = mean(collect(skipmissing(ARGO_OUTPUT[key]["Sₗ"][find_inverted])))
    pₘ = 0.5 .* (collect(skipmissing(ARGO_OUTPUT[key]["pₗ"][find_inverted])) .+
                collect(skipmissing(ARGO_OUTPUT[key]["pᵤ"][find_inverted])))
    pₘ_mean = mean(pₘ)
    αₗ = gsw_alpha.(Sₗ_mean, Θₗ_range, pₘ_mean)
    βₗ = gsw_beta.(Sₗ_mean, Θₗ_range, pₘ_mean)
    Δρ_thres = @. gsw_rho.(Sₗ_mean - (αₗ / βₗ) * ΔΘ_thres[i], Θₗ_range - ΔΘ_thres[i],
                        pₘ_mean) - gsw_rho.(Sₗ_mean, Θₗ_range, pₘ_mean)
    Δρ_thres_mean = mean(Δρ_thres)
    hist_edges = minimum(Δρˢ):bin_width:maximum(Δρˢ)

    hist!(ax[i], Δρˢ; bins = hist_edges, normalization = :pdf, color = (colours[i], 0.5))
    vlines!(ax[i], Δρ_thres_mean; color = :black, linestyle = :dash,
            label = "Δρ threshold for ΔΘ")
    ax[i].title = "PDF for ΔΘ = $(ΔΘ_thres[i])"
    hist_fit = fit(Histogram, Δρˢ, hist_edges)
    hist_fit = normalize(hist_fit; mode = :pdf)

    find_thres = findall(hist_edges .≤ Δρ_thres_mean)

    # to average threshold
    less_thres[i] = sum(hist_fit.weights[find_thres] .* bin_width)
    # after average threshold
    over_thres[i] = 1 - less_thres[i]
end
linkxaxes!(ax[1], ax[2])
linkxaxes!(ax[3], ax[4])
axislegend(ax[1]; position = :lt)
fig
less_thres

## density or ecdf plot that shows all four on same figure
fig2 = Figure(size = (500, 500))
ax = Axis(fig2[1, 1];
          title = "Empirical cumulative distribution for all ΔΘ's",
          xlabel = "Δρ (kgm⁻³)")
for (i, key) ∈ enumerate(ΔΘ_keys)
    Θₗ = collect(skipmissing(ARGO_OUTPUT[key]["Θₗ"]))
    Θᵤ = collect(skipmissing(ARGO_OUTPUT[key]["Θᵤ"]))
    find_inverted = findall(Θᵤ .< Θₗ)
    Δρˢ = collect(skipmissing(ARGO_OUTPUT[key]["Δρˢ"][find_inverted]))

    # density!(ax, Δρˢ; normalization = :pdf, color = (colours[i], 0.3),
    #          strokecolor = colours[i], strokearound = true, strokewidth = 3,
    #          label = "ΔΘ = $(ΔΘ_thres[i])")
   ecdfplot!(ax, Δρˢ; label = "ΔΘ = $(ΔΘ_thres[i])", color = colours[i])
end
vlines!(ax, 0; label = "Static instability", color = :black, linestyle = :dash)
axislegend(ax; position = :lt)
fig2

##
close(ARGO_OUTPUT)
