using .VerticalProfileStability
using JLD2, Statistics

extracted_data = jldopen(joinpath(@__DIR__, "ECCO_invertedΔΘ_extracted_data.jld2"))

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
colors = [:blue, :orange, :red, :green]
for (i, key) ∈ enumerate(keys(extracted_data))

    Θₗ, Δρˢ = extracted_data[key]["Θₗ"], extracted_data[key]["Δρˢ"]
    find = findall(xlimits[1] .≤ Θₗ .≤ xlimits[2] .&& ylimits[1] .≤ Δρˢ .≤ ylimits[2])
    Θₗ, Δρˢ = Θₗ[find], Δρˢ[find]
    lats = extracted_data[key]["lats"][find]
    Δρ_thres = extracted_data[key]["Δρ_thres"]
    ΔΘ_range = extracted_data[key]["ΔΘ_range"]
    Θ_lower_range = extracted_data[key]["Θ_lower_range"]

    sc = scatter!(ax, Θₗ, Δρˢ; color = colors[i], markersize = 4)
    lines!(ax, Θ_lower_range, Δρ_thres; color = colors[i],
           label = "ΔΘ = $(ΔΘ_range[1])ᵒC")

end
#vlines!(ax, Θ_lower_range; label = "Histogram bins", linestyle = :dash, color = :black)
#axislegend(ax; position = :rb)
Legend(fig[2, 1], ax, "Δρ threshold for", orientation = :horizontal)
save(joinpath(plotdir, "ECCO", "Θ_inversion", "2007_ΔΘ_thres_all.png"), fig)

## Full plot
# Temperature colourbar is wrong here, need a rethink for this and the pressure difference
keys_mat = reshape(keys(extracted_data), 2, 2)
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

    Θₗ, Δρˢ = extracted_data[key]["Θₗ"], extracted_data[key]["Δρˢ"]
    find = findall(xlimits[1] .≤ Θₗ .≤ xlimits[2] .&& ylimits[1] .≤ Δρˢ .≤ ylimits[2])
    Θₗ, Δρˢ = Θₗ[find], Δρˢ[find]
    plot_var = round.(extracted_data[key][colourbar_vars[choose_var]][find]; digits = 1)
    cbar_lab =  if colourbar_vars[choose_var] == "lats"
                   "Latitude (°N)"
                elseif colourbar_vars[choose_var] == "Δp_vals"
                    "Δp (dbar)"
                elseif colourbar_vars[choose_var] == "ΔΘ_vals"
                    "ΔΘ (°C)"
                end
    Δρ_thres = extracted_data[key]["Δρ_thres"]
    ΔΘ_range = extracted_data[key]["ΔΘ_range"]
    Θ_lower_range = extracted_data[key]["Θ_lower_range"]

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
save(joinpath(plotdir, "ECCO", "Θ_inversion",
              "2007_mult_ΔΘ_thres_$(colourbar_vars[choose_var]).png"), fig)

## Close data file
close(extracted_data)

## Statistics
extracted_data = jldopen(joinpath(@__DIR__, "ECCO_extracted_data.jld2"))
ΔΘ_0_5, ΔΘ_1, ΔΘ_2, ΔΘ_3 = keys(extracted_data)

## ΔΘ = 2 threshold
Θₗ = collect(skipmissing(extracted_data[ΔΘ_1]["Θₗ"]))
Δρˢ = collect(skipmissing(extracted_data[ΔΘ_1]["Δρˢ"]))
lats = extracted_data[ΔΘ_1]["lats"]
Δρ_thres = extracted_data[ΔΘ_1]["Δρ_thres"]
Θ_lower_range = extracted_data[ΔΘ_1]["Θ_lower_range"]

## Histogram
find_ = findall(Θₗ .≤ 10)
Θₗ_data = Θₗ[find_]
Δρˢ_data = Δρˢ[find_]
hist(Θₗ_data; bins = Θ_lower_range, normalization = :pdf)
hist(Θₗ_data; bins = range(extrema(Θ_lower_range)...; step = 0.5))

using StatsBase

## Fit a `Histogram` then get the index of the elements in each bin. This will allow
# comparison of `Δρ_thres` value calculated with all the static density differences for a
# given lower layer temperature.

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

##############################
## Individual plots if needed
for (i, key) ∈ enumerate(keys(extracted_data))

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

    Θₗ, Δρˢ = extracted_data[key]["Θₗ"], extracted_data[key]["Δρˢ"]
    find = findall(xlimits[1] .≤ Θₗ .≤ xlimits[2] .&& ylimits[1] .≤ Δρˢ .≤ ylimits[2])
    Θₗ, Δρˢ = Θₗ[find], Δρˢ[find]
    lats = extracted_data[key]["lats"][find]
    #ΔΘ = round.(extracted_data[key]["ΔΘ_vals"][find]; digits = 1)
    #Δp = round.(Int, extracted_data[key]["Δp_vals"][find])
    Δρ_thres = extracted_data[key]["Δρ_thres"]
    ΔΘ_range = extracted_data[key]["ΔΘ_range"]
    Θ_lower_range = extracted_data[key]["Θ_lower_range"]

    sc = scatter!(ax, Θₗ, Δρˢ; color = lats, markersize = 4)
    lines!(ax, Θ_lower_range, Δρ_thres; color = :red,
           label = "Density difference threshold for ΔΘ = $(ΔΘ_range[1])°C")

    Colorbar(fig[2, 1], sc, label = "Latitiude (°N)",
             vertical = false, flipaxis = false)

    axislegend(ax; position = :rb)

    @info "Saving file"
    save(joinpath(plotdir, "ECCO", "Θ_inversion", "2007_ΔΘ_thres_$(ΔΘ_range).png"), fig)
end
