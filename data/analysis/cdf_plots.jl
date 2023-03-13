# cdf plots and probabilities calculated for the ECCO and ARGO datasets.
using .VerticalProfileStability
using Statistics, StatsBase, ColorSchemes

## data files
const EXTRACTED_DATA_INV = joinpath(@__DIR__, "ECCO_invertedΔΘ_extracted_data.jld2")
const ARGO_OUTPUT = joinpath(@__DIR__, "ARGO_extracted.jld2")

data_files = (EXTRACTED_DATA_INV, ARGO_OUTPUT)
ΔΘ_thres = (0.5, 1, 2, 3)
ΔΘ_colours = get(ColorSchemes.thermal, range(0, 0.8; length = 4))

## setup figure
fig = Figure(size = (1000, 600))
titles = ("(a) CDF for ECCO data", "(b) CDF for ARGO data")
ax = [Axis(fig[i, 1];
           title = titles[i],
           xlabel = "Δρ (kgm⁻³)",
           ylabel = "ℙ(Δρₘˢ < Δρ | ΔΘ)") for i ∈ eachindex(data_files)]
xlims!(ax[1], -0.1, 0.01)
linkxaxes!(ax[1], ax[2])

Δρ_val = -0.04
ECCO_cdf_Δρ_val = Vector{Float64}(undef, 4)
ARGO_cdf_Δρ_val = Vector{Float64}(undef, 4)
## Calculate and plot ecdf
for (i, data) ∈ enumerate(data_files)

    d = jldopen(data)
    for (j, key) ∈ enumerate(keys(d))

        Δρˢ = sort(collect(skipmissing(d[key]["Δρˢ"])))
        fit_ecdf = ecdf(Δρˢ)
        lines!(ax[i], Δρˢ, fit_ecdf(Δρˢ);
              color = ΔΘ_colours[j],
              label = "ΔΘ = $(ΔΘ_thres[j])°C")
        i==1 ? ECCO_cdf_Δρ_val[j] = fit_ecdf(Δρ_val) :
               ARGO_cdf_Δρ_val[j] = fit_ecdf(Δρ_val)

    end
    vlines!(ax[i], Δρ_val; linestyle = :dash, color = :black, label = "Δρ' = $(Δρ_val)")
    close(d)

end

xlims!(ax[1], -0.5, 0.02)
xlims!(ax[2], -0.5, 0.02)
Legend(fig[3, :], ax[1], orientation = :horizontal)
fig
save(joinpath(PLOTDIR, "cdf.png"), fig)
ECCO_cdf_Δρ_val
ARGO_cdf_Δρ_val
##
test = jldopen(data_files[2])
Δρˢ = collect(skipmissing(test["ΔΘ_thres_[0.5, 1.0]"]["Δρˢ"]))
test_ecdf = ecdf(Δρˢ)
lines(sort(Δρˢ), test_ecdf(sort(Δρˢ)); color = :orange)
test_ecdf(-0.02)
close(test)