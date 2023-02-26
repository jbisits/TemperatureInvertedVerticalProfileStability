using .VerticalProfileStability

en4_data = readdir(EN4_DATADIR)[2]
en4_data_2021 = glob("*.nc", joinpath(EN4_DATADIR, en4_data))

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
