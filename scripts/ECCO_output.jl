## Saving files of variables from the monthly data
const ECCO_data_analysis = joinpath(@__DIR__, "../data/analysis/")
using Glob, Dates, JLD2
using .VerticalProfileStability

ECCO_data = glob("*.nc", ECCO_datadir)
timestamps = Date(2007, 01, 01):Day(1):Date(2007, 12, 31)
series = RasterSeries(ECCO_data, Ti(timestamps); child = RasterStack)

## ΔΘ_thres = 0.5
ΔΘ_thres = 0.5
mkdir(joinpath(ECCO_data_analysis, "output_$ΔΘ_thres"))
savepath = joinpath(ECCO_data_analysis, "output_$ΔΘ_thres")
series_max_Δρ(series, ΔΘ_thres, savepath)

## ΔΘ_thres = 1.0
ΔΘ_thres = 1.0
mkdir(joinpath(ECCO_data_analysis, "output_$ΔΘ_thres"))
savepath = joinpath(ECCO_data_analysis, "output_$ΔΘ_thres")
series_max_Δρ(series, ΔΘ_thres, savepath)

## ΔΘ_thres = 2.0
ΔΘ_thres = 2.0
mkdir(joinpath(ECCO_data_analysis, "output_$ΔΘ_thres"))
savepath = joinpath(ECCO_data_analysis, "output_$ΔΘ_thres")
series_max_Δρ(series, ΔΘ_thres, savepath)

## ΔΘ_thres = 3.0
ΔΘ_thres = 3.0
mkdir(joinpath(ECCO_data_analysis, "output_$ΔΘ_thres"))
savepath = joinpath(ECCO_data_analysis, "output_$ΔΘ_thres")
series_max_Δρ(series, ΔΘ_thres, savepath)

## 0.5 ≤ ΔΘ_thres < 1.0
ΔΘ_thres = [0.5, 1.0]
mkdir(joinpath(ECCO_data_analysis, "output_$ΔΘ_thres"))
savepath = joinpath(ECCO_data_analysis, "output_$ΔΘ_thres")
series_max_Δρ(series, ΔΘ_thres, savepath)

## 1.0 ≤ ΔΘ_thres < 2.0
ΔΘ_thres = [1.0, 2.0]
mkdir(joinpath(ECCO_data_analysis, "output_$ΔΘ_thres"))
savepath = joinpath(ECCO_data_analysis, "output_$ΔΘ_thres")
series_max_Δρ(series, ΔΘ_thres, savepath)

## 2.0 ≤ ΔΘ_thres < 3.0
ΔΘ_thres = [2.0, 3.0]
mkdir(joinpath(ECCO_data_analysis, "output_$ΔΘ_thres"))
savepath = joinpath(ECCO_data_analysis, "output_$ΔΘ_thres")
series_max_Δρ(series, ΔΘ_thres, savepath)

## Extract information from each data threshold and save to .jld2 for easy access
timestamps = Date(2007, 01, 01):Day(1):Date(2007, 12, 31)
ΔΘ_thres = [[0.5, 1.0], [1.0, 2.0], [2.0, 3.0]]
extracted_data = joinpath(ECCO_data_analysis, "ECCO_extracted_data.jld2")
for select_ΔΘ ∈ ΔΘ_thres
    output_path = joinpath(ECCO_data_analysis, "output_$(select_ΔΘ)")
    output_files = glob("*.nc", output_path)
    output_series = RasterSeries(output_files, Ti(timestamps); child = RasterStack)

    # Lower level Θ
    Θₗ = series2vec(output_series, :Θₗ)
    # Upper level Θ
    Θᵤ = series2vec(output_series, :Θᵤ)
    find_inversion = findall(Θᵤ .< Θₗ)
    # Extract all variables needed where there is a temperature inverison
    Θₗ = Θₗ[find_inversion]
    Θᵤ = Θᵤ[find_inversion]
    Sₗ = series2vec(output_series, :Sₗ)[find_inversion]
    pₗ = series2vec(output_series, :pₗ)[find_inversion]
    lats = get_lats(output_series, :Θₗ)[find_inversion]
    # Temperature and pressure differences
    ΔΘ_vals = series2vec(output_series, :ΔΘ)[find_inversion]
    Δp_vals = series2vec(output_series, :Δp)[find_inversion]
    # Density differneces
    Δρˢ = series2vec(output_series, :Δρ_static)[find_inversion]
    Δρᶜ = series2vec(output_series, :Δρ_cab)[find_inversion]

    # Density difference threshold
    Sₗ_mean = mean(Sₗ)
    pₗ_mean = mean(pₗ)

    Θ_lower_range = range(-1.85, 10; length = 100)
    Sₗ_mean_vec = fill(Sₗ_mean, 100)
    pₗ_mean_vec = fill(pₗ_mean, 100)
    α_vec = gsw_alpha.(Sₗ_mean_vec, Θ_lower_range , pₗ_mean_vec)
    β_vec = gsw_beta.(Sₗ_mean_vec, Θ_lower_range , pₗ_mean_vec)
    slope = α_vec ./ β_vec
    Δρ_thres =  gsw_rho.(Sₗ_mean_vec .- slope .* select_ΔΘ[1],
                        Θ_lower_range .- select_ΔΘ[1], pₗ_mean_vec) -
                gsw_rho.(Sₗ_mean_vec, Θ_lower_range, pₗ_mean_vec)

    jldopen(extracted_data, "a+") do file

        file["ΔΘ_thres_$(select_ΔΘ)"] = Dict("Θₗ" => Θₗ, "Δρˢ" => Δρˢ, "Δρᶜ" => Δρᶜ,
                                             "Δρ_thres" => Δρ_thres, "lats" => lats,
                                             "ΔΘ_range" => select_ΔΘ, "ΔΘ_vals" => ΔΘ_vals,
                                             "Δp_vals" => Δp_vals)

    end

end
