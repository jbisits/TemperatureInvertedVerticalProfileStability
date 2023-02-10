using .VerticalProfileStability
using JLD2, Statistics

extracted_data = jldopen(joinpath(@__DIR__, "ECCO_extracted_data.jld2"))

## Plots

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
    xlims!(ax, -1.88, 10)
    ylims!(ax, -0.1, 0.01)

    Θₗ, Δρˢ = extracted_data[key]["Θₗ"], extracted_data[key]["Δρˢ"]
    lats = extracted_data[key]["lats"]
    Δρ_thres = extracted_data[key]["Δρ_thres"]
    ΔΘ_range = extracted_data[key]["ΔΘ_range"]
    Θ_lower_range = range(-1.85, 10; length = 100) # forgot to save this

    sc = scatter!(ax, Θₗ, Δρˢ; color = lats, markersize = 4)
    lines!(ax, Θ_lower_range, Δρ_thres; color = :red,
           label = "Density difference threshold for ΔΘ = $(ΔΘ_range)(ᵒC)")

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
xlims!(ax, -1.88, 10)
ylims!(ax, -0.1, 0.01)
colors = [:blue, :orange, :red]
for (i, key) ∈ enumerate(keys(extracted_data))

    Θₗ, Δρˢ = extracted_data[key]["Θₗ"], extracted_data[key]["Δρˢ"]
    lats = extracted_data[key]["lats"]
    Δρ_thres = extracted_data[key]["Δρ_thres"]
    ΔΘ_range = extracted_data[key]["ΔΘ_range"]
    Θ_lower_range = range(-1.85, 10; length = 100) # forgot to save this

    sc = scatter!(ax, Θₗ, Δρˢ; color = colors[i], markersize = 4)
    lines!(ax, Θ_lower_range, Δρ_thres; color = colors[i],
           label = "Density difference threshold for ΔΘ = $(ΔΘ_range[1])ᵒC")

end
axislegend(ax; position = :rb)
save(joinpath(plotdir, "ECCO", "2007_ΔΘ_thres_all.png"), fig)
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
    xlims!(ax[i], -1.88, 10)
    ylims!(ax[i], -0.1, 0.01)
end
fig

## Plot
for (i, key) ∈ enumerate(keys(extracted_data))

    Θₗ, Δρˢ = extracted_data[key]["Θₗ"], extracted_data[key]["Δρˢ"]
    lats = extracted_data[key]["lats"]
    Δρ_thres = extracted_data[key]["Δρ_thres"]
    ΔΘ_range = extracted_data[key]["ΔΘ_range"]
    Θ_lower_range = range(-1.85, 10; length = 100) # forgot to save this

    sc = scatter!(ax, Θₗ, Δρˢ; color = lats, markersize = 4)
    lines!(ax, Θ_lower_range, Δρ_thres; color = :red,
           label = "Density difference threshold for ΔΘ = $(ΔΘ_range)(ᵒC)")

    if i == length(keys(extracted_data))
        Colorbar(fig[1, 4], sc, label = "Latitude (ᵒN)", vertical = false, flipaxis = false)
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
Θ_lower_range = range(-1.85, 10; length = 100) # forgot to save this
close(extracted_data)

## Loess, dont think this is that useful
find_ = findall(Θₗ .< 10 .&& Δρˢ .> -0.04)
perm = sortperm(Θₗ[find_])
xs = Float64.(Θₗ[find_][perm])
ys = Float64.(Δρˢ[find_][perm])
using Loess

model = loess(xs, ys; span = 0.25)
us = range(extrema(xs)...; length = 200)
vs = predict(model, us)
lines(us, vs)
