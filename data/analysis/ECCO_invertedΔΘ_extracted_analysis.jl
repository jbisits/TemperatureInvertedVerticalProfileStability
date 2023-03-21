using .VerticalProfileStability
using JLD2, Statistics, ColorSchemes, StatsBase, LinearAlgebra

const EXTRACTED_DATA_INV = joinpath(@__DIR__, "ECCO_invertedΔΘ_extracted_data.jld2")
inv_data = jldopen(EXTRACTED_DATA_INV)

## Plots
# set ylimits, much faster to extract then plot data then plot and use lims! on axis.
xlimits = (-1.88, 10)
ylimits = (-0.1, 0.01)

## Same axis, thresholds in different colours
fig = Figure(size = (700, 700))
ax = Axis(fig[1, 1];
        xlabel = "Θ (ᵒC)",
        xaxisposition = :top,
        title = "Maximum static density difference between two vertically\nspaced levels of a profile against temperature of lower level",
        ylabel = "Δρ (kgm⁻³)")
xlims!(ax, xlimits)
ylims!(ax, ylimits)
#colours = [:blue, :orange, :red, :green]
colours = get(ColorSchemes.thermal, range(0, 0.8; length = 4))
for (i, key) ∈ enumerate(keys(inv_data))

    Θₗ, Δρˢ = inv_data[key]["Θₗ"], inv_data[key]["Δρˢ"]
    find = findall(xlimits[1] .≤ Θₗ .≤ xlimits[2] .&& ylimits[1] .≤ Δρˢ .≤ ylimits[2])
    Θₗ, Δρˢ = Θₗ[find], Δρˢ[find]
    lats = inv_data[key]["lats"][find]
    Δρ_thres = inv_data[key]["Δρ_thres"]
    ΔΘ_range = inv_data[key]["ΔΘ_range"]
    Θ_lower_range = inv_data[key]["Θ_lower_range"]

    sc = scatter!(ax, Θₗ, Δρˢ; color = colours[i], markersize = 4)
    lines!(ax, Θ_lower_range, Δρ_thres; color = colours[i],
           label = "ΔΘ = $(ΔΘ_range[1])ᵒC", linewidth = 2)

end
#vlines!(ax, Θ_lower_range; label = "Histogram bins", linestyle = :dash, color = :black)
#axislegend(ax; position = :rb)
Legend(fig[2, 1], ax, "Δρ threshold for", orientation = :horizontal, nbanks = 2)
save(joinpath(PLOTDIR, "ECCO", "Θ_inversion", "2007_ΔΘ_thres_all.png"), fig)

## Full plot
# Temperature colourbar is wrong here, need a rethink for this and the pressure difference
keys_mat = reshape(keys(inv_data), 2, 2)
colourbar_vars = ("lats", "Δp_vals", "ΔΘ_vals")
colourbar_col = (:viridis, :batlow, :thermal)
# These have been found by looking at the subset for each key in the dictionary and taking
# the `extrema`. The code is not here though.
colourbar_ranges = ((-90, 90), (20, 410), (-3.5, -0.5))
choose_var = 3
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
    xlims!(ax[i, j], xlimits)
    ylims!(ax[i, j], ylimits)
end
fig

## Plot
for (i, key) ∈ enumerate(keys_mat)

    Θₗ, Δρˢ = inv_data[key]["Θₗ"], inv_data[key]["Δρˢ"]
    find = findall(xlimits[1] .≤ Θₗ .≤ xlimits[2] .&& ylimits[1] .≤ Δρˢ .≤ ylimits[2])
    Θₗ, Δρˢ = Θₗ[find], Δρˢ[find]
    plot_var = round.(inv_data[key][colourbar_vars[choose_var]][find]; digits = 1)
    cbar_lab =  if colourbar_vars[choose_var] == "lats"
                   "Latitude (°N)"
                elseif colourbar_vars[choose_var] == "Δp_vals"
                    "Δp (dbar)"
                elseif colourbar_vars[choose_var] == "ΔΘ_vals"
                    "ΔΘ (°C)"
                end
    Δρ_thres = inv_data[key]["Δρ_thres"]
    ΔΘ_range = inv_data[key]["ΔΘ_range"]
    Θ_lower_range = inv_data[key]["Θ_lower_range"]

    sc = scatter!(ax[i], Θₗ, Δρˢ;
                  markersize = 4,
                  color = plot_var,
                  colormap = colourbar_col[choose_var],
                  colorrange = colourbar_ranges[choose_var])
    lines!(ax[i], Θ_lower_range, Δρ_thres;
           color = colors[i],
           label = "ΔΘ = $(ΔΘ_range[1])ᵒC")

    if i == 4
        Colorbar(fig[:, 3], sc, label = cbar_lab)
    end
end
lin_elements = [LineElement(color = col) for col ∈ colors]
leg_labels = "ΔΘ = " .* string.([0.5, 1.0, 2.0, 3.0]) .* "°C"
Legend(fig[3, :], lin_elements, leg_labels, "Δρ threshold for", orientation = :horizontal)
# Save
save(joinpath(PLOTDIR, "ECCO", "Θ_inversion",
              "2007_mult_ΔΘ_thres_$(colourbar_vars[choose_var]).png"), fig)

## Close data file
close(inv_data)
##

ΔΘ_keys = keys(inv_data)
num_obs = Array{Int64}(undef, length(ΔΘ_keys))
for (i, key) ∈ enumerate(ΔΘ_keys)
    num_obs[i] = length(inv_data[key]["Θₗ"])
end
num_obs
bin_width = 0.01

fig = Figure(size = (1200, 1200))
ax = [Axis(fig[i, j];
        xlabel = "Δρ (kgm⁻³)"
        )
      for i ∈ 1:2, j ∈ 1:2]
less_thres = Vector{Float64}(undef, 4)
over_thres = Vector{Float64}(undef, 4)

for (i, key) ∈ enumerate(ΔΘ_keys)

    Θₗ, Δρˢ = inv_data[key]["Θₗ"], inv_data[key]["Δρˢ"]
    Δρ_thres = inv_data[key]["Δρ_thres"]

    Δρ_thres_mean = mean(Δρ_thres)
    hist_edges = minimum(Δρˢ):bin_width:maximum(Δρˢ)

    hist!(ax[i], Δρˢ; bins = hist_edges, normalization = :pdf, color = (colours[i], 0.5))
    vlines!(ax[i], Δρ_thres_mean; color = colours[i], linestyle = :dash,
            label = "Δρ threshold for ΔΘ", linewidth = 2)
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
# zoom in
for i ∈ 1:4
    xlims!(ax[i], -1, 0.01)
end
fig
## density or ecdf plot that shows all four on same figure
fig2 = Figure(size = (500, 500))
ax = Axis(fig2[1, 1];
          title = "Empirical cumulative distribution for all ΔΘ's",
          xlabel = "Δρ (kgm⁻³)")
for (i, key) ∈ enumerate(ΔΘ_keys)

    Θₗ, Δρˢ = collect(skipmissing(inv_data[key]["Θₗ"])),
              collect(skipmissing(inv_data[key]["Δρˢ"]))

    Δρ_thres = inv_data[key]["Δρ_thres"]
    Δρ_thres_mean = mean(Δρ_thres)
    # density!(ax, Δρˢ; normalization = :pdf, color = (colours[i], 0.3),
    #          strokecolor = colours[i], strokearound = true, strokewidth = 3,
    #          label = "ΔΘ = $(ΔΘ_thres[i])")
   ecdfplot!(ax, Δρˢ; label = "ΔΘ = $(ΔΘ_thres[i])", color = colours[i])
   vlines!(ax, Δρ_thres_mean; color = colours[i], linestyle = :dash)
end
vlines!(ax, 0; label = "Static instability", color = :black, linestyle = :dash)
axislegend(ax; position = :lt)
fig2

## zoom in
xlims!(ax, -1, 0.01)
fig2
############################################################################################
## full fig - scatter and pdfs.
############################################################################################
Δρ_lims = (-0.1, 0.01)
xlimits = (-1.88, 6)
ylimits = (-0.1, 0.01)
bin_width = 0.0001
full_fig = Figure(size = (800, 1000))
colours = get(ColorSchemes.thermal, range(0, 0.8; length = 4))
# scatter
splot = full_fig[1, 1] = GridLayout()
ax_splot = Axis(splot[1, 1];
                xlabel = "Θ (ᵒC)",
                xaxisposition = :top,
                title = "(a) Maximum static density difference between two vertically\nspaced levels of a profile against temperature of lower level",
                ylabel = "Δρ (kgm⁻³)")
for (i, key) ∈ enumerate(keys(inv_data))

    Θₗ, Δρˢ = inv_data[key]["Θₗ"], inv_data[key]["Δρˢ"]
    lats = inv_data[key]["lats"]
    find = findall(xlimits[1] .≤ Θₗ .≤ xlimits[2] .&& ylimits[1] .≤ Δρˢ .≤ ylimits[2] .&&
                   lats .≤ -60)
    Θₗ, Δρˢ = Θₗ[find], Δρˢ[find]
    Θ_lower_range = inv_data[key]["Θ_lower_range"]
    find_below_6 = findall(Θ_lower_range .≤ 6)
    Θ_lower_range = inv_data[key]["Θ_lower_range"][find_below_6]
    Δρ_thres = inv_data[key]["Δρ_thres"][find_below_6]
    ΔΘ_range = inv_data[key]["ΔΘ_range"]

    sc = scatter!(ax_splot, Θₗ, Δρˢ; color = colours[i], markersize = 4)
    lines!(ax_splot, Θ_lower_range, Δρ_thres; color = colours[i],
           label = "ΔΘ = $(ΔΘ_range[1])ᵒC", linewidth = 2)

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

for (i, key) ∈ enumerate(keys(inv_data))

    Θₗ, Δρˢ = inv_data[key]["Θₗ"], inv_data[key]["Δρˢ"]
    lats = inv_data[key]["lats"]
    find = findall(lats .≤ -60)
    Θₗ, Δρˢ = Θₗ[find], Δρˢ[find]
    ΔΘ_range = inv_data[key]["ΔΘ_range"]
    Δρ_thres = inv_data[key]["Δρ_thres"]

    Δρ_thres_mean = mean(Δρ_thres)
    hist_edges = minimum(Δρˢ):bin_width:maximum(Δρˢ)

    hist!(ax_pdf[i], Δρˢ; bins = hist_edges, normalization = :pdf, color = (colours[i], 0.5))
    vlines!(ax_pdf[i], Δρ_thres_mean; color = colours[i], linewidth = 2,
            label = "Δρ threshold for ΔΘ")
    vlines!(ax_pdf[i], 0; color = :black, linestyle = :dash)
    ax_pdf[i].title = letter_labels[i] * " PDF for ΔΘ = $(ΔΘ_range[1])°C"
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
#full_fig
save(joinpath(PLOTDIR, "ECCO", "Θ_inversion", "ecco_sc_pdf_60S.png"), full_fig)
##
close(inv_data)
## Below here maybe not that useful.

## Statistics
ΔΘ_0_5, ΔΘ_1, ΔΘ_2, ΔΘ_3 = keys(inv_data)

## ΔΘ = 2 threshold
Θₗ = collect(skipmissing(inv_data[ΔΘ_3]["Θₗ"]))
Δρˢ = collect(skipmissing(inv_data[ΔΘ_3]["Δρˢ"]))
lats = inv_data[ΔΘ_3]["lats"]
Δρ_thres = inv_data[ΔΘ_3]["Δρ_thres"]
mean(Δρ_thres)
Θ_lower_range = inv_data[ΔΘ_3]["Θ_lower_range"]

## Histogram
find_ = findall(Θₗ .≤ 10)
Θₗ_data = Θₗ[find_]
Δρˢ_data = Δρˢ[find_]
hist(Θₗ_data; bins = Θ_lower_range, normalization = :pdf)
hist(Θₗ_data; bins = range(extrema(Θ_lower_range)...; step = 0.5))

## Fit a `Histogram` then get the index of the elements in each bin. This will allow
# comparison of `Δρ_thres` value calculated with all the static density differences for a
# given lower layer temperature.
using StatsBase

Θₗ_histfit = fit(Histogram, Θₗ_data, Θ_lower_range)
Θₗ_binidx = StatsBase.binindex.(Ref(Θₗ_histfit), Θₗ_data)
unique_binidx = unique(Θₗ_binidx)
Δρ_thres_binidx = Δρ_thres[unique_binidx]
unstable_pts = Vector{Float64}(undef, length(unique_binidx))
for (i, idx) ∈ enumerate(unique_binidx)

    find = findall(idx .== Θₗ_binidx)
    unstable_pts[i] = length(findall(Δρˢ_data[find] .> Δρ_thres_binidx[i]))

end
unstable_pts
sum(unstable_pts)
find_unstable_idx = findall(unstable_pts .> 0)
unstable_Θₗ = Θ_lower_range[find_unstable_idx]
sum(unstable_pts) / length(Δρˢ_data)

Θ_lower_range[unique_binidx[67]]
find = findall(unique_binidx[67] .== Θₗ_binidx)
length(findall(Δρˢ_data[find] .> Δρ_thres[67]))

############################################################################################
## Individual plots if needed
for (i, key) ∈ enumerate(keys(inv_data))

    @info "Generating plot for $(key)"
    fig = Figure(size = (600, 600))
    ax = Axis(fig[1, 1];
            xlabel = "Θ (ᵒC)",
            xaxisposition = :top,
            title = "Maximum density difference between two vertically spaced\nlevels of a profile against temperature of lower level",
            ylabel = "Δρ (kgm⁻³)",
            subtitle = key*"ᵒC")
    xlims!(ax, xlimits)
    ylims!(ax, ylimits)

    Θₗ, Δρˢ = inv_data[key]["Θₗ"], inv_data[key]["Δρˢ"]
    find = findall(xlimits[1] .≤ Θₗ .≤ xlimits[2] .&& ylimits[1] .≤ Δρˢ .≤ ylimits[2])
    Θₗ, Δρˢ = Θₗ[find], Δρˢ[find]
    lats = inv_data[key]["lats"][find]
    #ΔΘ = round.(inv_data[key]["ΔΘ_vals"][find]; digits = 1)
    #Δp = round.(Int, inv_data[key]["Δp_vals"][find])
    Δρ_thres = inv_data[key]["Δρ_thres"]
    ΔΘ_range = inv_data[key]["ΔΘ_range"]
    Θ_lower_range = inv_data[key]["Θ_lower_range"]

    sc = scatter!(ax, Θₗ, Δρˢ; color = lats, markersize = 4)
    lines!(ax, Θ_lower_range, Δρ_thres; color = :red,
           label = "Density difference threshold for ΔΘ = $(ΔΘ_range[1])°C")

    Colorbar(fig[2, 1], sc, label = "Latitiude (°N)",
             vertical = false, flipaxis = false)

    axislegend(ax; position = :rb)

    @info "Saving file"
    save(joinpath(PLOTDIR, "ECCO", "Θ_inversion", "2007_ΔΘ_thres_$(ΔΘ_range).png"), fig)
end
close(inv_data)
