## Saving files of variables from the monthly data
const ECCO_DATA_ANALYSIS = joinpath(@__DIR__, "../data/analysis/")
using Glob, Dates, JLD2
using .VerticalProfileStability

ECCO_data = glob("*.nc", ECCO_DATADIR)
timestamps = Date(2007, 01, 01):Day(1):Date(2007, 12, 31)
rs_series = RasterSeries(ECCO_data, Ti(timestamps); child = RasterStack)

## ΔΘ_thres = 0.5
ΔΘ_thres = 0.5
mkdir(joinpath(ECCO_DATA_ANALYSIS, "output_$ΔΘ_thres"))
savepath = joinpath(ECCO_DATA_ANALYSIS, "output_$ΔΘ_thres")
series_max_Δρ(rs_series, ΔΘ_thres, savepath)

## ΔΘ_thres = 1.0
ΔΘ_thres = 1.0
mkdir(joinpath(ECCO_DATA_ANALYSIS, "output_$ΔΘ_thres"))
savepath = joinpath(ECCO_DATA_ANALYSIS, "output_$ΔΘ_thres")
series_max_Δρ(rs_series, ΔΘ_thres, savepath)

## ΔΘ_thres = 2.0
ΔΘ_thres = 2.0
mkdir(joinpath(ECCO_DATA_ANALYSIS, "output_$ΔΘ_thres"))
savepath = joinpath(ECCO_DATA_ANALYSIS, "output_$ΔΘ_thres")
series_max_Δρ(rs_series, ΔΘ_thres, savepath)

## ΔΘ_thres = 3.0
ΔΘ_thres = 3.0
mkdir(joinpath(ECCO_DATA_ANALYSIS, "output_$ΔΘ_thres"))
savepath = joinpath(ECCO_DATA_ANALYSIS, "output_$ΔΘ_thres")
series_max_Δρ(rs_series, ΔΘ_thres, savepath)

## 0.5 ≤ ΔΘ_thres < 1.0
ΔΘ_thres = [0.5, 1.0]
mkdir(joinpath(ECCO_DATA_ANALYSIS, "output_$ΔΘ_thres"))
savepath = joinpath(ECCO_DATA_ANALYSIS, "output_$ΔΘ_thres")
series_max_Δρ(rs_series, ΔΘ_thres, savepath)

## 1.0 ≤ ΔΘ_thres < 2.0
ΔΘ_thres = [1.0, 2.0]
mkdir(joinpath(ECCO_DATA_ANALYSIS, "output_$ΔΘ_thres"))
savepath = joinpath(ECCO_DATA_ANALYSIS, "output_$ΔΘ_thres")
series_max_Δρ(rs_series, ΔΘ_thres, savepath)

## 2.0 ≤ ΔΘ_thres < 3.0
ΔΘ_thres = [2.0, 3.0]
mkdir(joinpath(ECCO_DATA_ANALYSIS, "output_$ΔΘ_thres"))
savepath = joinpath(ECCO_DATA_ANALYSIS, "output_$ΔΘ_thres")
series_max_Δρ(rs_series, ΔΘ_thres, savepath)

## Extract temperature inverted data from each data threshold and save to .jld2
timestamps = Date(2007, 01, 01):Day(1):Date(2007, 12, 31)
ΔΘ_thres = [[0.5, 1.0], [1.0, 2.0], [2.0, 3.0], 3.0]
extracted_data = joinpath(ECCO_DATA_ANALYSIS, "ECCO_invertedΔΘ_extracted_data.jld2")
for select_ΔΘ ∈ ΔΘ_thres
    # Data
    @info "Reading RasterSeries"
    output_path = joinpath(ECCO_DATA_ANALYSIS, "output_$(select_ΔΘ)")
    output_files = glob("*.nc", output_path)
    output_series = RasterSeries(output_files, Ti(timestamps); child = RasterStack)

    @info "Extacting temperature inverted data"
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
    pᵤ = series2vec(output_series, :pᵤ)[find_inversion]
    p̄ = @. 0.5 * (pᵤ + pₗ)
    lats = get_lats(output_series, :Θₗ)[find_inversion]

    # Temperature and pressure differences
    ΔΘ_vals = series2vec(output_series, :ΔΘ)[find_inversion]
    Δp_vals = series2vec(output_series, :Δp)[find_inversion]

    # Density differneces
    Δρˢ = series2vec(output_series, :Δρ_static)[find_inversion]
    Δρᶜ = series2vec(output_series, :Δρ_cab)[find_inversion]

    # Density difference threshold
    @info "Computing density difference threshold"
    Sₗ_mean = mean(Sₗ)
    p̄_mean = mean(p̄)
    Θ_lower_range = range(-1.85, 10; length = 100)
    Sₗ_mean_vec = fill(Sₗ_mean, 100)
    p̄_mean_vec = fill(p̄_mean, 100)
    α_vec = gsw_alpha.(Sₗ_mean_vec, Θ_lower_range , p̄_mean_vec)
    β_vec = gsw_beta.(Sₗ_mean_vec, Θ_lower_range , p̄_mean_vec)
    slope = α_vec ./ β_vec
    ΔΘ = typeof(select_ΔΘ) == Vector{Float64} ? select_ΔΘ[1] : select_ΔΘ
    Δρ_thres =  gsw_rho.(Sₗ_mean_vec .- slope .* ΔΘ, Θ_lower_range .- ΔΘ, p̄_mean_vec) -
                gsw_rho.(Sₗ_mean_vec, Θ_lower_range, p̄_mean_vec)

    @info "Saving $(select_ΔΘ)"
    jldopen(extracted_data, "a+") do file

        file["ΔΘ_thres_$(select_ΔΘ)"] = Dict("Θₗ" => Θₗ, "Δρˢ" => Δρˢ, "Δρᶜ" => Δρᶜ,
                                             "Δρ_thres" => Δρ_thres, "lats" => lats,
                                             "Θ_lower_range" => Θ_lower_range,
                                             "ΔΘ_range" => select_ΔΘ, "ΔΘ_vals" => ΔΘ_vals,
                                             "Δp_vals" => Δp_vals)

    end

end

## Extract temperature stratified data from each data threshold and save to .jld2
timestamps = Date(2007, 01, 01):Day(1):Date(2007, 12, 31)
ΔΘ_thres = [[0.5, 1.0], [1.0, 2.0], [2.0, 3.0]]
extracted_data = joinpath(ECCO_DATA_ANALYSIS, "ECCO_stratifiedΔΘ_extracted_data.jld2")
for select_ΔΘ ∈ ΔΘ_thres
    # Data
    @info "Reading RasterSeries"
    output_path = joinpath(ECCO_DATA_ANALYSIS, "output_$(select_ΔΘ)")
    output_files = glob("*.nc", output_path)
    output_series = RasterSeries(output_files, Ti(timestamps); child = RasterStack)

    @info "Extacting temperature stratified data"
    # Lower level Θ
    Θₗ = series2vec(output_series, :Θₗ)
    # Upper level Θ
    Θᵤ = series2vec(output_series, :Θᵤ)
    find_stratified = findall(Θᵤ .> Θₗ)

    # Extract all variables needed where there is a temperature inverison
    Θₗ = Θₗ[find_stratified]
    Θᵤ = Θᵤ[find_stratified]
    Sₗ = series2vec(output_series, :Sₗ)[find_stratified]
    pₗ = series2vec(output_series, :pₗ)[find_stratified]
    pᵤ = series2vec(output_series, :pᵤ)[find_stratified]
    p̄ = @. 0.5 * (pᵤ + pₗ)
    lats = get_lats(output_series, :Θₗ)[find_stratified]

    # Temperature and pressure differences
    ΔΘ_vals = series2vec(output_series, :ΔΘ)[find_stratified]
    Δp_vals = series2vec(output_series, :Δp)[find_stratified]

    # Density differneces
    Δρˢ = series2vec(output_series, :Δρ_static)[find_stratified]
    Δρᶜ = series2vec(output_series, :Δρ_cab)[find_stratified]

    # Density difference threshold
    @info "Computing density difference threshold"
    Sₗ_mean = mean(Sₗ)
    p̄_mean = mean(p̄)
    Θ_lower_range = range(-1.85, 10; length = 100)
    Sₗ_mean_vec = fill(Sₗ_mean, 100)
    p̄_mean_vec = fill(p̄_mean, 100)
    α_vec = gsw_alpha.(Sₗ_mean_vec, Θ_lower_range , p̄_mean_vec)
    β_vec = gsw_beta.(Sₗ_mean_vec, Θ_lower_range , p̄_mean_vec)
    slope = α_vec ./ β_vec
    ΔΘ = typeof(select_ΔΘ) == Vector{Float64} ? select_ΔΘ[1] : select_ΔΘ
    Δρ_thres =  gsw_rho.(Sₗ_mean_vec .+ slope .* ΔΘ,
                         Θ_lower_range .+ ΔΘ, p̄_mean_vec) -
                gsw_rho.(Sₗ_mean_vec, Θ_lower_range, p̄_mean_vec)

    @info "Saving $(select_ΔΘ)"
    jldopen(extracted_data, "a+") do file

        file["ΔΘ_thres_$(select_ΔΘ)"] = Dict("Θₗ" => Θₗ, "Δρˢ" => Δρˢ, "Δρᶜ" => Δρᶜ,
                                             "Δρ_thres" => Δρ_thres, "lats" => lats,
                                             "Θ_lower_range" => Θ_lower_range,
                                             "ΔΘ_range" => select_ΔΘ, "ΔΘ_vals" => ΔΘ_vals,
                                             "Δp_vals" => Δp_vals)

    end

end

## Temperature inverted profile location map

output_path = joinpath(ECCO_DATA_ANALYSIS, "output_[0.5, 1.0]")
output_files = glob("*.nc", output_path)
ΔΘ_series = RasterSeries(output_files, Dim{:dates}(timestamps);
                         child = Raster, name = "ΔΘ")

for rs ∈ ΔΘ_series
    set(rs, )
    map!(x -> ismissing(x) ? x : x ≤ -0.5 ? x : missing, rs, rs)
end

heatmap(ΔΘ_series[end])
length(ΔΘ_series)

## Change to 1, then can merge and sum over the time dimension.

for rs ∈ ΔΘ_series
    map!(x -> ismissing(x) ? 0 : x ≤ -0.5 ? 1 : 0, rs, rs)
end

ΔΘ_counts =  ΔΘ_series[1].data[:, :, 1]
for rs ∈ ΔΘ_series[2:end]
    ΔΘ_counts .+= rs.data[:, :, 1]
end
map!(x -> x == 0 ? missing : x, ΔΘ_counts, ΔΘ_counts)

ΔΘ_proportion = ΔΘ_counts ./ sum(skipmissing(ΔΘ_counts))
heatmap(ΔΘ_proportion)

rs_proportion = Raster(ΔΘ_proportion, (X(lookup(ΔΘ_series[1], :X)), Y(lookup(ΔΘ_series[1], :Y))))

using GeoMakie
fig = Figure(size = (500, 500))
ax = GeoAxis(fig[1, 1];
             title = "Location and count of temperature\ninverted profiles for ECCO",
             xlabel = "Longitude",
             ylabel = "Latitude",
             coastlines = true)
hm = heatmap!(ax, lookup(rs_proportion, :X), lookup(rs_proportion, :Y), ΔΘ_counts;
              colormap = :batlow)
#Colorbar(fig[1, 2], hm, label = "Number of profiles")
Colorbar(fig[2, 1], hm, label = "Number of profiles", vertical = false, flipaxis = true)
fig
#save(joinpath(PLOTDIR, "ECCO/Θ_inversion/profile_count.png"), fig)
# create full figure by running corresponding geo makie map in goship_output.jl
#save(joinpath(PLOTDIR, "combined_profile_location.png"), fig)
