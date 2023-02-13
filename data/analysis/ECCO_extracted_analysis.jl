using .VerticalProfileStability
using JLD2, Statistics

extracted_data = jldopen(joinpath(@__DIR__, "ECCO_extracted_data.jld2"))

## Plots
# set ylimits, much faster to extract then plot data then plot and use lims! on axis.
xlimits = (-1.88, 10)
ylimits = (-0.1, 0.01)
## Plot, individual plots as full plot took over an hour before I gave up waiting
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
    Δρ_thres = extracted_data[key]["Δρ_thres"]
    ΔΘ_range = extracted_data[key]["ΔΘ_range"]
    Θ_lower_range = extracted_data[key]["Θ_lower_range"]

    sc = scatter!(ax, Θₗ, Δρˢ; color = lats, markersize = 4)
    lines!(ax, Θ_lower_range, Δρ_thres; color = :red,
           label = "Density difference threshold for ΔΘ = $(ΔΘ_range[1])(ᵒC)")

    Colorbar(fig[2, 1], sc, label = "Latitude (ᵒN)", vertical = false, flipaxis = false)

    axislegend(ax; position = :rb)

    @info "Saving file"
    save(joinpath(plotdir, "ECCO", "2007_ΔΘ_thres_$(ΔΘ_range).png"), fig)
end

## Same axis, thresholds in different colours
fig = Figure(size = (600, 600))
ax = Axis(fig[1, 1];
        xlabel = "Θ (ᵒC)",
        xaxisposition = :top,
        title = "Maximum density difference between two vertically spaced\nlevels of a profile against temperature of lower level",
        ylabel = "Δρ (kgm⁻³)")
xlims!(ax, xlimits)
ylims!(ax, ylimits)
colors = [:blue, :orange, :red]
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
           label = "Density difference threshold for ΔΘ = $(ΔΘ_range[1])ᵒC")

end
vlines!(ax, Θ_lower_range; label = "Histogram bins", linestyle = :dash, color = :black)
axislegend(ax; position = :rb)
save(joinpath(plotdir, "ECCO", "2007_ΔΘ_thres_all_withbins.png"), fig)
## Full plot

## Setup axis
fig = Figure(size = (600, 1400))
ax = [Axis(fig[i, 1];
          xlabel = "Θ (ᵒC)",
          xaxisposition = :top,
          ylabel = "Δρ (kgm⁻³)",
          subtitle = ΔΘ*"ᵒC")
      for (i, ΔΘ) ∈ enumerate(keys(extracted_data))]
ax[1].title = "Maximum density difference between two vertically spaced\nlevels of a profile against temperature of lower level"
linkxaxes!(ax...)
hidexdecorations!(ax[2], grid = false)
hidexdecorations!(ax[3], grid = false)
for i ∈ eachindex(keys(extracted_data))
    xlims!(ax[i], xlimits)
    ylims!(ax[i], ylimits)
end
fig

## Plot
for (i, key) ∈ enumerate(keys(extracted_data))

    Θₗ, Δρˢ = extracted_data[key]["Θₗ"], extracted_data[key]["Δρˢ"]
    find = findall(xlimits[1] .≤ Θₗ .≤ xlimits[2] .&& ylimits[1] .≤ Δρˢ .≤ ylimits[2])
    Θₗ, Δρˢ = Θₗ[find], Δρˢ[find]
    lats = extracted_data[key]["lats"][find]
    Δρ_thres = extracted_data[key]["Δρ_thres"]
    ΔΘ_range = extracted_data[key]["ΔΘ_range"]
    Θ_lower_range = extracted_data[key]["Θ_lower_range"]

    sc = scatter!(ax[i], Θₗ, Δρˢ; color = lats, markersize = 4)
    lines!(ax[i], Θ_lower_range, Δρ_thres; color = :red,
           label = "Density difference threshold for ΔΘ = $(ΔΘ_range)(ᵒC)")

    if i == length(keys(extracted_data))
        Colorbar(fig[4, 1], sc, label = "Latitude (ᵒN)", vertical = false, flipaxis = false)
    end

    axislegend(ax[i]; position = :rb)
end

## Save
save(joinpath(plotdir, "ECCO", "2007_mult_ΔΘ_thres.png"), fig)
## Close data file
close(extracted_data)

## Statistics

extracted_data = jldopen(joinpath(@__DIR__, "ECCO_extracted_data.jld2"))
data_keys = keys(extracted_data)

## ΔΘ = 2 threshold

Θₗ = collect(skipmissing(extracted_data[data_keys[3]]["Θₗ"]))
Δρˢ = collect(skipmissing(extracted_data[data_keys[3]]["Δρˢ"]))
lats = extracted_data[data_keys[3]]["lats"]
Δρ_thres = extracted_data[data_keys[3]]["Δρ_thres"]
Θ_lower_range = extracted_data[data_keys[3]]["Θ_lower_range"]
close(extracted_data)

## Find proportion above theoretical threshold
find_ = findall(Θₗ .≤ 10)
Θₗ_data = Θₗ[find_]
Δρˢ_data = Δρˢ[find_]
hist(Θₗ_data; bins = Θ_lower_range)

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

Θ_lower_range[unique_binidx[68]]
find = findall(unique_binidx[67] .== Θₗ_binidx)
length(findall(Δρˢ_data[find] .> Δρ_thres[67]))


## Loess, dont think this is that useful
# find_ = findall(Θₗ .< 10 .&& Δρˢ .> -0.04)
# perm = sortperm(Θₗ[find_])
# xs = Float64.(Θₗ[find_][perm])
# ys = Float64.(Δρˢ[find_][perm])
# using Loess

# model = loess(xs, ys; span = 0.25)
# us = range(extrema(xs)...; length = 200)
# vs = predict(model, us)
# lines(us, vs)
