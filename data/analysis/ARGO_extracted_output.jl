using .VerticalProfileStability
using Statistics

const ARGO_OUTPUT = joinpath(@__DIR__, "ARGO_extracted.jld2")
argo_data = jldopen(ARGO_OUTPUT)

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
#colours = [:blue, :orange, :red, :green]
colours = get(ColorSchemes.thermal, range(0, 0.8; length = 4))
for (i, key) ∈ enumerate(keys(argo_data))

    Θᵤ = collect(skipmissing(argo_data[key]["Θᵤ"]))
    Θₗ = collect(skipmissing(argo_data[key]["Θₗ"]))
    find_inverted = findall(Θᵤ .< Θₗ)
    Δρˢ = collect(skipmissing(argo_data[key]["Δρˢ"]))
    sc = scatter!(ax, Θₗ[find_inverted], Δρˢ[find_inverted];
                  color = colours[i], markersize = 5)

    Sₗ_mean = mean(collect(skipmissing(argo_data[key]["Sₗ"][find_inverted])))
    pₘ = 0.5 .* (collect(skipmissing(argo_data[key]["pₗ"][find_inverted])) .+
                 collect(skipmissing(argo_data[key]["pᵤ"][find_inverted])))
    pₘ_mean = mean(pₘ)
    αₗ = gsw_alpha.(Sₗ_mean, Θₗ_range, pₘ_mean)
    βₗ = gsw_beta.(Sₗ_mean, Θₗ_range, pₘ_mean)
    Δρ_thres = @. gsw_rho.(Sₗ_mean - (αₗ / βₗ) * ΔΘ_thres[i], Θₗ_range - ΔΘ_thres[i],
                          pₘ_mean) - gsw_rho.(Sₗ_mean, Θₗ_range, pₘ_mean)
    lines!(ax, Θₗ_range, Δρ_thres; color = colours[i],
           label = "ΔΘ = $(ΔΘ_thres[i])ᵒC")

end
#axislegend(ax; position = :rb)
Legend(fig[2, 1], ax, "Δρ threshold for", orientation = :horizontal)
fig
#save(joinpath(PLOTDIR, "ARGO", "ΔΘ_thres_all.png"), fig)
ylims!(ax, -0.1, 0.01)
fig
save(joinpath(PLOTDIR, "ARGO", "ΔΘ_thres_all_zoom.png"), fig)
## Full plot
# Temperature colourbar is wrong here, need a rethink for this and the pressure difference
keys_mat = reshape(keys(argo_data), 2, 2)
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

    Θᵤ = collect(skipmissing(argo_data[key]["Θᵤ"]))
    Θₗ = collect(skipmissing(argo_data[key]["Θₗ"]))
    Δρˢ = collect(skipmissing(argo_data[key]["Δρˢ"]))
    plot_var = round.(argo_data[key][colourbar_vars[choose_var]][find]; digits = 1)
    find_inverted = findall(Θᵤ .< Θₗ)
    # plot_var = round.(argo_data[key][colourbar_vars[choose_var]][find_inverted]; digits = 1)
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

ΔΘ_keys = keys(argo_data)
bin_width = 0.01

fig = Figure(size = (1200, 1200))
ax = [Axis(fig[i, j];
        xlabel = "Δρ (kgm⁻³)"
        )
      for i ∈ 1:2, j ∈ 1:2]
less_thres = Vector{Float64}(undef, 4)
over_thres = Vector{Float64}(undef, 4)

for (i, key) ∈ enumerate(ΔΘ_keys)
    Θₗ = collect(skipmissing(argo_data[key]["Θₗ"]))
    Θᵤ = collect(skipmissing(argo_data[key]["Θᵤ"]))
    find_inverted = findall(Θᵤ .< Θₗ)
    Δρˢ = collect(skipmissing(argo_data[key]["Δρˢ"][find_inverted]))
    Sₗ_mean = mean(collect(skipmissing(argo_data[key]["Sₗ"][find_inverted])))
    pₘ = 0.5 .* (collect(skipmissing(argo_data[key]["pₗ"][find_inverted])) .+
                collect(skipmissing(argo_data[key]["pᵤ"][find_inverted])))
    pₘ_mean = mean(pₘ)
    αₗ = gsw_alpha.(Sₗ_mean, Θₗ_range, pₘ_mean)
    βₗ = gsw_beta.(Sₗ_mean, Θₗ_range, pₘ_mean)
    Δρ_thres = @. gsw_rho.(Sₗ_mean - (αₗ / βₗ) * ΔΘ_thres[i], Θₗ_range - ΔΘ_thres[i],
                        pₘ_mean) - gsw_rho.(Sₗ_mean, Θₗ_range, pₘ_mean)
    Δρ_thres_mean = mean(Δρ_thres)
    hist_edges = minimum(Δρˢ):bin_width:maximum(Δρˢ)

    hist!(ax[i], Δρˢ; bins = hist_edges, normalization = :pdf, color = (colours[i], 0.5))
    vlines!(ax[i], Δρ_thres_mean; color = colours[i], linestyle = :dash,
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
for i ∈ 1:4
    xlims!(ax[i], -0.2, 0.01)
end
fig
## density or ecdf plot that shows all four on same figure
fig2 = Figure(size = (500, 500))
ax = Axis(fig2[1, 1];
          title = "Empirical cumulative distribution for all ΔΘ's",
          xlabel = "Δρ (kgm⁻³)")
for (i, key) ∈ enumerate(ΔΘ_keys)
    Θₗ = collect(skipmissing(argo_data[key]["Θₗ"]))
    Θᵤ = collect(skipmissing(argo_data[key]["Θᵤ"]))
    find_inverted = findall(Θᵤ .< Θₗ)
    Δρˢ = collect(skipmissing(argo_data[key]["Δρˢ"][find_inverted]))

    # density!(ax, Δρˢ; normalization = :pdf, color = (colours[i], 0.3),
    #          strokecolor = colours[i], strokearound = true, strokewidth = 3,
    #          label = "ΔΘ = $(ΔΘ_thres[i])")
   ecdfplot!(ax, Δρˢ; label = "ΔΘ = $(ΔΘ_thres[i])", color = colours[i])
end
vlines!(ax, 0; label = "Static instability", color = :black, linestyle = :dash)
axislegend(ax; position = :lt)
fig2

############################################################################################
## full fig - scatter and pdfs.
############################################################################################
Δρ_lims = (-0.2, 0.01)
bin_width = 0.0025
full_fig = Figure(size = (800, 1000))
# scatter
splot = full_fig[1, 1] = GridLayout()
ax_splot = Axis(splot[1, 1];
                xlabel = "Θ (ᵒC)",
                xaxisposition = :top,
                title = "(a) Maximum static density difference between two vertically\nspaced levels of a profile against temperature of lower level",
                ylabel = "Δρ (kgm⁻³)")
for (i, key) ∈ enumerate(keys(argo_data))

    Θᵤ = collect(skipmissing(argo_data[key]["Θᵤ"]))
    Θₗ = collect(skipmissing(argo_data[key]["Θₗ"]))
    find_inverted = findall(Θᵤ .< Θₗ)
    Δρˢ = collect(skipmissing(argo_data[key]["Δρˢ"]))
    sc = scatter!(ax_splot, Θₗ[find_inverted], Δρˢ[find_inverted];
                    color = colours[i], markersize = 5)

    Sₗ_mean = mean(collect(skipmissing(argo_data[key]["Sₗ"][find_inverted])))
    pₘ = 0.5 .* (collect(skipmissing(argo_data[key]["pₗ"][find_inverted])) .+
                    collect(skipmissing(argo_data[key]["pᵤ"][find_inverted])))
    pₘ_mean = mean(pₘ)
    αₗ = gsw_alpha.(Sₗ_mean, Θₗ_range, pₘ_mean)
    βₗ = gsw_beta.(Sₗ_mean, Θₗ_range, pₘ_mean)
    Δρ_thres = @. gsw_rho.(Sₗ_mean - (αₗ / βₗ) * ΔΘ_thres[i], Θₗ_range - ΔΘ_thres[i],
                            pₘ_mean) - gsw_rho.(Sₗ_mean, Θₗ_range, pₘ_mean)
    lines!(ax_splot, Θₗ_range, Δρ_thres; color = colours[i], linewidth = 2,
            label = "ΔΘ = $(ΔΘ_thres[i])ᵒC")

end
ylims!(ax_splot, Δρ_lims)
Legend(splot[1, 2], ax_splot, "Δρ threshold for")

# pdf
pdf_plots = full_fig[2:3, 1] = GridLayout()
ax_pdf = [Axis(pdf_plots[i, j];
        xlabel = "Δρ (kgm⁻³)"
        )
      for i ∈ 1:2, j ∈ 1:2]
less_thres = Vector{Float64}(undef, 4)
over_thres = Vector{Float64}(undef, 4)
letter_labels = ["(b)", "(c)", "(d)", "(e)"]

for (i, key) ∈ enumerate(keys(argo_data))
    Θₗ = collect(skipmissing(argo_data[key]["Θₗ"]))
    Θᵤ = collect(skipmissing(argo_data[key]["Θᵤ"]))
    find_inverted = findall(Θᵤ .< Θₗ)
    Δρˢ = collect(skipmissing(argo_data[key]["Δρˢ"][find_inverted]))
    Sₗ_mean = mean(collect(skipmissing(argo_data[key]["Sₗ"][find_inverted])))
    pₘ = 0.5 .* (collect(skipmissing(argo_data[key]["pₗ"][find_inverted])) .+
                collect(skipmissing(argo_data[key]["pᵤ"][find_inverted])))
    pₘ_mean = mean(pₘ)
    αₗ = gsw_alpha.(Sₗ_mean, Θₗ_range, pₘ_mean)
    βₗ = gsw_beta.(Sₗ_mean, Θₗ_range, pₘ_mean)
    Δρ_thres = @. gsw_rho.(Sₗ_mean - (αₗ / βₗ) * ΔΘ_thres[i], Θₗ_range - ΔΘ_thres[i],
                        pₘ_mean) - gsw_rho.(Sₗ_mean, Θₗ_range, pₘ_mean)
    Δρ_thres_mean = mean(Δρ_thres)
    hist_edges = minimum(Δρˢ):bin_width:maximum(Δρˢ)

    hist!(ax_pdf[i], Δρˢ; bins = hist_edges, normalization = :pdf, color = (colours[i], 0.5))
    vlines!(ax_pdf[i], Δρ_thres_mean; color = colours[i], linewidth = 2,
            label = "Δρ threshold for ΔΘ")
    vlines!(ax_pdf[i], 0; color = :black, linestyle = :dash,
            label = "Static stability threshold")
    ax_pdf[i].title = letter_labels[i] * " PDF for ΔΘ = $(ΔΘ_thres[i])°C"
    hist_fit = fit(Histogram, Δρˢ, hist_edges)
    hist_fit = normalize(hist_fit; mode = :pdf)

    find_thres = findall(hist_edges .≤ Δρ_thres_mean)

    # to average threshold
    less_thres[i] = sum(hist_fit.weights[find_thres] .* bin_width)
    # after average threshold
    over_thres[i] = 1 - less_thres[i]

end
pdf_lims = (-0.2, 0.01)
for i ∈ 1:4
    xlims!(ax_pdf[i], pdf_lims)
end
less_thres
over_thres
rowsize!(full_fig.layout, 1, Auto(1.15))
full_fig

##
save(joinpath(PLOTDIR, "ARGO", "argo_ss_pdf.png"), full_fig)
##
close(argo_data)
