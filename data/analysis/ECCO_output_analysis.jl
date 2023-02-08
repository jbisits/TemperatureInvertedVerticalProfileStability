using .VerticalProfileStability
using Statistics

timestamps = Date(2007, 01, 01):Day(1):Date(2007, 12, 31)
ΔΘ_thres = [0.5, 1.0, 2.0, [1.0, 2.5]]
select_ΔΘ = 4
output_path = joinpath(@__DIR__, "output_$(ΔΘ_thres[select_ΔΘ])")
output_files = glob("*.nc", output_path)
output_series = RasterSeries(output_files, Ti(timestamps); child = RasterStack)

#########
## Note
# I had the order wrong when saving to symbols. This means that
# :ΔΘ => :pᵤ, :pᵤ => :pₗ, :pₗ => :ΔΘ.
# This has been corrected in the function now but means might need to re run
# so that code can be reused. No just change the symbols in the `RasterSeries` for
# the ΔΘ = 2 output.
########

Θₗ = series2vec(output_series, :Θₗ)
lats = get_lats(output_series, :Θₗ)
Δρˢ = series2vec(output_series, :Δρ_static)
Δρᶜ = series2vec(output_series, :Δρ_cab)
Δρ = [Δρˢ, Δρᶜ]

## Plot, only look at Δρ = [-0.1, 0.01] and Θ = (-1.89, 5) to start
fig = Figure(size = (1200, 600))
dd = ["static", "cabbeling"]
ax = [Axis(fig[1, i];
          title = "Maximum $(dd[i]) density difference of a\nprofile against temperature of lower level\nwhere max dd was calculated with ΔΘ = $(ΔΘ_thres[select_ΔΘ])",
          subtitle = "ECCO data 2007",
          xlabel = "Θ (ᵒC) at lower level",
          xaxisposition = :top,
          ylabel = "Maximum Δρ (kgm⁻³) of profile") for i ∈ 1:2]
linkyaxes!(ax...)
hideydecorations!(ax[2], grid = false)
for i ∈ 1:2
    xlims!(ax[i], -1.88, 10)
    ylims!(ax[i], -0.1, maximum(Δρᶜ))
end
for (i, Δρ_) ∈ enumerate(Δρ)
    sc = scatter!(ax[i], Θₗ, Δρ_; markersize = 4, color = lats)
    #scatter!(ax[i], Θₗ, Δρ; markersize = 4) # no colour
    lines!(ax[i], zeros(4), 0:-1:-3;
           linestyle = :dash, linewidth = 3, color = :red)
    lines!(ax[i], 0.5 .* ones(4), 0:-1:-3;
           linestyle = :dash, linewidth = 3, color = :orange)
    lines!(ax[i], -1.88 .* ones(4), 0:-1:-3;
           linestyle = :dash, linewidth = 3, color = :green)
    if i == length(Δρ)
        Colorbar(fig[2, :], sc, label = "Latitude (ᵒN)", vertical = false, flipaxis = false)
    end
end
#fig
save(joinpath(plotdir, "ECCO", "2007_ΔΘ_thres_$(ΔΘ_thres[select_ΔΘ])_zoom.png"), fig)

## Restrict to temperature inversion
Θᵤ = series2vec(output_series, :Θᵤ)
find_inversion = findall(Θᵤ .< Θₗ )

Θₗ_inversion = Θₗ[find_inversion]
Δρ_inversion = [Δρˢ[find_inversion], Δρᶜ[find_inversion]]

fig = Figure(size = (1200, 600))
dd = ["static", "cabbeling"]
ax = [Axis(fig[1, i];
          title = "Maximum $(dd[i]) density difference of a\nprofile against temperature of lower level\nwhere max dd was calculated with ΔΘ = $(ΔΘ_thres[select_ΔΘ])",
          subtitle = "ECCO data 2007",
          xlabel = "Θ (ᵒC) at lower level",
          xaxisposition = :top,
          ylabel = "Maximum Δρ (kgm⁻³) of profile") for i ∈ 1:2]
linkyaxes!(ax...)
hideydecorations!(ax[2], grid = false)
for i ∈ 1:2
    xlims!(ax[i], -1.88, 10)
    ylims!(ax[i], -0.1, maximum(Δρᶜ))
end
for (i, Δρ_) ∈ enumerate(Δρ_inversion)
    sc = scatter!(ax[i], Θₗ_inversion, Δρ_; markersize = 4, color = lats[find_inversion])
    #scatter!(ax[i], Θₗ, Δρ; markersize = 4) # no colour
    lines!(ax[i], zeros(4), 0:-1:-3;
           linestyle = :dash, linewidth = 3, color = :red)
    lines!(ax[i], 0.5 .* ones(4), 0:-1:-3;
           linestyle = :dash, linewidth = 3, color = :orange)
    lines!(ax[i], -1.88 .* ones(4), 0:-1:-3;
           linestyle = :dash, linewidth = 3, color = :green)
    if i == length(Δρ)
        Colorbar(fig[2, :], sc, label = "Latitude (ᵒN)", vertical = false, flipaxis = false)
    end
end
#fig

## Theoretical curve
Sₗ_mean = mean(series2vec(output_series, :Sₗ))
pₗ_mean = mean(series2vec(output_series, :pₗ))

# Loook at temperature values
ΔΘ_vals = series2vec(output_series, :ΔΘ)
minimum(ΔΘ_vals), maximum(ΔΘ_vals)
Θₗ_inversion_mean = mean(Θₗ_inversion)
Θᵤ_inversion_mean = mean(Θᵤ[find_inversion])
ΔΘ_inversion_mean = mean(ΔΘ_vals[find_inversion])
abs(Θᵤ_inversion_mean - Θₗ_inversion_mean)

Θ_lower_range = range(-1.85, 10; length = 100)
Sₗ_mean_vec = fill(Sₗ_mean, 100)
pₗ_mean_vec = fill(pₗ_mean, 100)
α_vec = gsw_alpha.(Sₗ_mean_vec, Θ_lower_range , pₗ_mean_vec)
β_vec = gsw_beta.(Sₗ_mean_vec, Θ_lower_range , pₗ_mean_vec)
slope = α_vec ./ β_vec
Δρ_thres_u =  gsw_rho.(Sₗ_mean_vec .- slope, Θ_lower_range .- 1, pₗ_mean_vec) -
              gsw_rho.(Sₗ_mean_vec, Θ_lower_range, pₗ_mean_vec)
Δρ_thres_l =  gsw_rho.(Sₗ_mean_vec .- slope .* 2, Θ_lower_range .- 2, pₗ_mean_vec) -
              gsw_rho.(Sₗ_mean_vec, Θ_lower_range, pₗ_mean_vec)
lines!(ax[1], Θ_lower_range, Δρ_thres_u; color = :red)
lines!(ax[1], Θ_lower_range, Δρ_thres_l; color = :red)
fig
