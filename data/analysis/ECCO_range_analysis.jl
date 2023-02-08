using .VerticalProfileStability
using Statistics

## Read in the data and extract vectors
timestamps = Date(2007, 01, 01):Day(1):Date(2007, 12, 31)
ΔΘ_thres = [[0.5, 1.0], [1.0, 2.0], [2.0, 3.0]]
select_ΔΘ = 3
output_path = joinpath(@__DIR__, "output_$(ΔΘ_thres[select_ΔΘ])")
output_files = glob("*.nc", output_path)
output_series = RasterSeries(output_files, Ti(timestamps); child = RasterStack)

Θₗ = series2vec(output_series, :Θₗ)
Θᵤ = series2vec(output_series, :Θᵤ)
find_inversion = findall(Θᵤ .< Θₗ )
Sₗ = series2vec(output_series, :Sₗ)
pₗ = series2vec(output_series, :pₗ)
lats = get_lats(output_series, :Θₗ)
Δρˢ = series2vec(output_series, :Δρ_static)
Δρᶜ = series2vec(output_series, :Δρ_cab)
Δρ = [Δρˢ, Δρᶜ]

Θₗ_inversion = Θₗ[find_inversion]
Δρ_inversion = [Δρˢ[find_inversion], Δρᶜ[find_inversion]]


## Look at temperature inverted profiles and values
ΔΘ_vals = series2vec(output_series, :ΔΘ)
(minimum(ΔΘ_vals), maximum(ΔΘ_vals))
Θₗ_inversion_mean = mean(Θₗ_inversion)
Θᵤ_inversion_mean = mean(Θᵤ[find_inversion])
ΔΘ_inversion_mean = mean(ΔΘ_vals[find_inversion])
abs(Θᵤ_inversion_mean - Θₗ_inversion_mean)

## Scatter plot
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

## Theoretical Δρ thresholds

Sₗ_mean = mean(Sₗ[find_inversion])
pₗ_mean = mean(pₗ[find_inversion])

Θ_lower_range = range(-1.85, 10; length = 100)
Sₗ_mean_vec = fill(Sₗ_mean, 100)
pₗ_mean_vec = fill(pₗ_mean, 100)
α_vec = gsw_alpha.(Sₗ_mean_vec, Θ_lower_range , pₗ_mean_vec)
β_vec = gsw_beta.(Sₗ_mean_vec, Θ_lower_range , pₗ_mean_vec)
slope = α_vec ./ β_vec
Δρ_thres_u =  gsw_rho.(Sₗ_mean_vec .- slope .* ΔΘ_thres[select_ΔΘ][1],
                       Θ_lower_range .- ΔΘ_thres[select_ΔΘ][1], pₗ_mean_vec) -
              gsw_rho.(Sₗ_mean_vec, Θ_lower_range, pₗ_mean_vec)
Δρ_thres_l =  gsw_rho.(Sₗ_mean_vec .- slope .* ΔΘ_thres[select_ΔΘ][2],
                       Θ_lower_range .- ΔΘ_thres[select_ΔΘ][2], pₗ_mean_vec) -
              gsw_rho.(Sₗ_mean_vec, Θ_lower_range, pₗ_mean_vec)
lines!(ax[1], Θ_lower_range, Δρ_thres_u; color = :red)
lines!(ax[1], Θ_lower_range, Δρ_thres_l; color = :red)
fig
