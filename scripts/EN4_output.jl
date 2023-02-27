using .VerticalProfileStability
using Statistics, JLD2

## Choose a threshold for all years
ΔΘ_thres = [0.5, 1.0, 2.0, 3.0]
en4_yearly_data = readdir(EN4_DATADIR)[2:end]
saved_output = joinpath(@__DIR__, "..", "data", "analysis", "EN4_extracted.jld2")
for year ∈ en4_yearly_data

    @info year[end-3:end] #* "ΔΘ_thres = "
    en4_data = glob("*.nc", joinpath(EN4_DATADIR, year))
    output = MaximumDensityDifference.en4_max_Δρ(en4_data, ΔΘ_thres[2])
    jldopen(saved_output, "a+") do file
        file["EN4_"*year[end-3:end]] = output
    end

end


## Look at output and check things are working for single year of data
test_en4_Δρ_max = MaximumDensityDifference.en4_max_Δρ(en4_yearly_data[1], 3.0)
Θᵤ = test_en4_Δρ_max["Θᵤ"]
Sᵤ = test_en4_Δρ_max["Sᵤ"]
pᵤ = test_en4_Δρ_max["pᵤ"]
Θₗ = test_en4_Δρ_max["Θₗ"]
Sₗ = test_en4_Δρ_max["Sₗ"]
pₗ = test_en4_Δρ_max["pₗ"]
Δρˢ = test_en4_Δρ_max["Δρˢ"]
find_inversion = findall(Θᵤ .< Θₗ) ∩ findall(Θₗ .≤ 10)

fig = Figure(size = (500, 500))
ax = Axis(fig[1, 1];
          xlabel = "Θₗ (°C)",
          ylabel = "Δρ (kgm⁻³)")
scatter!(ax, Θₗ[find_inversion], Δρˢ[find_inversion])
Θ_range = range(extrema(Θₗ[find_inversion])...; length = 100)
Sₗ_mean = mean(Sₗ[find_inversion])
p̄_mean = mean(0.5 .* (pₗ[find_inversion] .+ pᵤ[find_inversion]))
αₗ = gsw_alpha.(Sₗ_mean, Θ_range, p̄_mean)
βₗ = gsw_beta.(Sₗ_mean, Θ_range, p̄_mean)

Δρ_thres = @. gsw_rho(Sₗ_mean - (αₗ / βₗ) * 3, Θ_range - 3, p̄_mean) - gsw_rho.(Sₗ_mean, Θ_range, p̄_mean)
lines!(ax, Θ_range, Δρ_thres)
fig
## Exploring data and looking at best way to extract and calculate density difference.
ds = NCDataset(en4_data_2021; aggdim = "N_PROF")
lat = ds["LATITUDE"][:]
lon = ds["LONGITUDE"][:]
θ = ds["POTM_CORRECTED"][:, :]
Sₚ = ds["PSAL_CORRECTED"][:, :]
z = -ds["DEPH_CORRECTED"][:, :]

for i ∈ eachindex(lat)
    println(i)
    profile_depths = collect(skipmissing(z[:, i]))
    find_1000 = findall(profile_depths .≥ -1000)
    if !isnothing(find_1000)

        z_profile = nomissing(z[find_1000, i], NaN)
        Sₚ_profile = nomissing(Sₚ[find_1000, i], NaN)
        θ_profile = nomissing(θ[find_1000, i], NaN)
        p = gsw_p_from_z.(z_profile, lat[i])
        Sₐ = gsw_sa_from_sp.(Sₚ_profile, p, lon[i], lat[i])
        Θ = gsw_ct_from_pt.(θ_profile, Sₐ)

        res = Δρ_max(Sₐ, Θ, p, ΔΘ_thres)

    end

end

test_find_d = findall(collect(skipmissing(z[:, 872])) .>= -1000)
test_find_872 = findall(.!ismissing.(Sₚ[test_find_d, 872])) ∩ findall(.!ismissing.(θ[test_find_d, 872])) ∩ findall(.!ismissing.(z[test_find_d, 872]))
findall(ismissing(Sₚ[test_find_872, 872]))
findall(ismissing(θ[test_find_872, 872]))
θ[test_find_872, 872]
collect(skipmissing(z[test_find_872, 872])) == z[test_find_872, 872]
test_find = findall(.!ismissing.(Sₚ[:, 872])) ∩ findall(.!ismissing.(θ[:, 872]))
Sₚ[test_find, 872]
findall(ismissing(collect(Sₚ[test_find, 872])))
findall(ismissing(θ[test_find, 873]))
collect(Sₚ[test_find, 873])

findall(.!ismissing.(Sₚ[test_find, 872]))==
findall(.!ismissing.(θ[test_find, 872]))
findall(.!ismissing.(z[test_find, 872]))
