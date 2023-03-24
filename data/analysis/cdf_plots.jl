# cdf plots and probabilities calculated for the ECCO and ARGO datasets.
using .VerticalProfileStability
using Statistics, StatsBase, ColorSchemes

## data files
const EXTRACTED_DATA_INV = joinpath(@__DIR__, "ECCO_invertedΔΘ_extracted_data.jld2")
const ARGO_OUTPUT = joinpath(@__DIR__, "ARGO_extracted.jld2")
const GOSHIP_JOINED = joinpath(@__DIR__, "goship_joined.jld2")

data_files = (EXTRACTED_DATA_INV, ARGO_OUTPUT)
ΔΘ_thres = (0.5, 1, 2, 3)
ΔΘ_colours = get(ColorSchemes.thermal, range(0, 0.8; length = 4))

## setup figure, ECCO and ARGO
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
ECCO_num_obs = Vector{Int64}(undef, 4)
ARGO_cdf_Δρ_val = Vector{Float64}(undef, 4)
ARGO_num_obs = Vector{Int64}(undef, 4)
## Calculate and plot ecdf
for (i, data) ∈ enumerate(data_files)

    d = jldopen(data)
    for (j, key) ∈ enumerate(keys(d))

        lats = collect(skipmissing(d[key]["lats"]))
        find = findall(lats .≤ -60)
        Δρˢ = sort(collect(skipmissing(d[key]["Δρˢ"][find])))
        # Δρˢ = sort(collect(skipmissing(d[key]["Δρˢ"])))
        i==1 ? ECCO_num_obs[j] = length(Δρˢ) :
               ARGO_num_obs[j] = length(Δρˢ)
        fit_ecdf = ecdf(Δρˢ)
        lines!(ax[i], Δρˢ, fit_ecdf(Δρˢ);
              color = ΔΘ_colours[j],
              label = "ΔΘ = $(ΔΘ_thres[j])°C")
        i==1 ? ECCO_cdf_Δρ_val[j] = fit_ecdf(Δρ_val) :
               ARGO_cdf_Δρ_val[j] = fit_ecdf(Δρ_val)

    end
    vlines!(ax[i], 0; linestyle = :dash, color = :black)
    vlines!(ax[i], Δρ_val; linestyle = :dash, color = :red, label = "Δρ' = $(Δρ_val)")
    close(d)

end

xlims!(ax[1], -0.5, 0.02)
xlims!(ax[2], -0.5, 0.02)
Legend(fig[3, :], ax[1], orientation = :horizontal)
fig
save(joinpath(PLOTDIR, "cdf_60S.png"), fig)
ECCO_cdf_Δρ_val
ARGO_cdf_Δρ_val
ECCO_num_obs
ARGO_num_obs

## setup figure, ECCO and GOSHIP
data_files = (EXTRACTED_DATA_INV, GOSHIP_JOINED)
ΔΘ_thres = (0.5, 1, 2, 3)
ΔΘ_colours = get(ColorSchemes.thermal, range(0, 0.8; length = 4))

fig = Figure(size = (1000, 600))
titles = ("(a) CDF for ECCO data", "(b) CDF for GOSHIP data")
ax = [Axis(fig[i, 1];
           title = titles[i],
           xlabel = "Δρ (kgm⁻³)",
           ylabel = "ℙ(Δρₘˢ < Δρ | ΔΘ)") for i ∈ eachindex(data_files)]
xlims!(ax[1], -0.1, 0.01)
linkxaxes!(ax[1], ax[2])

Δρ_val = -0.04
ECCO_cdf_Δρ_val = Vector{Float64}(undef, 4)
ECCO_num_obs = Vector{Int64}(undef, 4)
GOSHIP_cdf_Δρ_val = Vector{Float64}(undef, 4)
GOSHIP_num_obs = Vector{Int64}(undef, 4)
## Calculate and plot ecdf
for (i, data) ∈ enumerate(data_files)

    d = jldopen(data)
    for (j, key) ∈ enumerate(keys(d))

        Δρˢ =   if i == 1
                    sort(collect(skipmissing(d[key]["Δρˢ"])))
                else
                    Θᵤ = collect(skipmissing(d[key]["Θᵤ"]))
                    Θₗ = collect(skipmissing(d[key]["Θₗ"]))
                    find_inv = Θᵤ .≤ Θₗ
                    sort(collect(skipmissing(d[key]["Δρˢ"]))[find_inv])
                end
        i==1 ? ECCO_num_obs[j] = length(Δρˢ) :
               GOSHIP_num_obs[j] = length(Δρˢ)
        fit_ecdf = ecdf(Δρˢ)
        lines!(ax[i], Δρˢ, fit_ecdf(Δρˢ);
              color = ΔΘ_colours[j],
              label = "ΔΘ = $(ΔΘ_thres[j])°C")
        i==1 ? ECCO_cdf_Δρ_val[j] = fit_ecdf(Δρ_val) :
               GOSHIP_cdf_Δρ_val[j] = fit_ecdf(Δρ_val)

    end
    vlines!(ax[i], 0; linestyle = :dash, color = :black)
    vlines!(ax[i], Δρ_val; linestyle = :dash, color = :red, label = "Δρ' = $(Δρ_val)")
    close(d)

end

xlims!(ax[1], -0.5, 0.02)
xlims!(ax[2], -0.5, 0.02)
Legend(fig[3, :], ax[1], orientation = :horizontal)
fig
save(joinpath(PLOTDIR, "cdf_ECCO_GOSHIP.png"), fig)
ECCO_cdf_Δρ_val
GOSHIP_cdf_Δρ_val
ECCO_num_obs
GOSHIP_num_obs
##

##
test = jldopen(data_files[2])
Δρˢ = collect(skipmissing(test["ΔΘ_thres_3.0"]["Δρˢ"]))
lats = collect(skipmissing(test["ΔΘ_thres_[0.5, 1.0]"]["lats"]))
find = findall(lats .<= -60)
test_ecdf = ecdf(Δρˢ[find])
lines(sort(Δρˢ), test_ecdf(sort(Δρˢ)); color = :orange)
test_ecdf(-0.02)
close(test)
