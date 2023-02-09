using .VerticalProfileStability
using JLD2


extracted_data = jldopen(joinpath(@__DIR__, "ECCO_extracted_data.jld2"))

## Setup axis
fig = Figure(size = (600, 1800))
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

    sc = scatter!(ax[i], extracted_data[key]["Θₗ"], extracted_data[key]["Δρˢ"];
                  color = extracted_data[key]["lats"],)
    lines!(ax[i], extracted_data[key]["Θₗ"], extracted_data[key]["Δρˢ"]; color = :red,
           label = "Density difference threshold for ΔΘ = $(extracted_data[key]["ΔΘ_range"])(ᵒC)")

    if i == length(keys(extracted_data))
        Colorbar(fig[1, 4], sc, label = "Latitude (ᵒN)", vertical = false, flipaxis = false)
    end

    axislegend(ax[i]; position = :rb)
end

## Save
save(joinpath(plotdir, "ECCO", "2007_mult_ΔΘ_thres.png"), fig)
## Close data file
close(extracted_data)
