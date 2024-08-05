using .VerticalProfileStability, ColorSchemes

## data
ECCO_files = glob("*.nc", ECCO_TS_DATA_PATH)
series_path = "/Users/Joey/Documents/PhD data and code/VerticalProfileStability/scripts/../data/analysis/output_1.0"
timestamps = Date(2007, 01, 01):Day(1):Date(2007, 12, 31)
rs_s = RasterSeries(series_path, Ti(timestamps), child = RasterStack)

extracted_data = jldopen(EXTRACTED_DATA_INV)
          Δρˢ = extracted_data["ΔΘ_thres_1.0"]["Δρˢ"]
          Δρᶜ = extracted_data["ΔΘ_thres_1.0"]["Δρᶜ"]
computed_lats = extracted_data["ΔΘ_thres_1.0"]["lats"]
Θ_lower_range = extracted_data["ΔΘ_thres_1.0"]["Θ_lower_range"]
     Δρ_thres = extracted_data["ΔΘ_thres_1.0"]["Δρ_thres"]
close(extracted_data)
## Filter profiles
Δρ_full = []
Θᵤ_full = []
Θₗ_full = []
lons_full = []
lats_full = []
dates_full = []

for t ∈ eachindex(timestamps)

    Δρ = rs_s[t][:Δρ_static].data
    replace!(Δρ, missing => NaN)
    Δρ_cab = rs_s[t][:Δρ_cab].data
    replace!(Δρ_cab, missing => NaN)
    Θᵤ = rs_s[t][:Θᵤ].data
    replace!(Θᵤ, missing => NaN)
    Θₗ = rs_s[t][:Θₗ].data
    replace!(Θₗ, missing => NaN)

    find = findall(Θᵤ .< Θₗ)
    lon = lookup(rs_s[1], X)
    lat = lookup(rs_s[1], Y)
    lon_grid = lon .* ones(length(lat))'
    lat_grid = ones(length(lon)) .* lat'

    Δρ_invΘ = vec(Δρ[find])
    Δρ_invΘ_cab = vec(Δρ_cab[find])
    Θₗ_ = vec(Θₗ[find])
    Θᵤ_ = vec(Θᵤ[find])
    longitudes = lon_grid[find]
    latidiudes = lat_grid[find]
    date = fill(timestamps[t], length(find))

    push!(Δρ_full, Δρ_invΘ...)
    push!(Θᵤ_full, Θᵤ_...)
    push!(Θₗ_full, Θₗ_...)
    push!(lons_full, longitudes...)
    push!(lats_full, latidiudes...)
    push!(dates_full, date...)

end

find_unstable = findall(Δρ_full .> 0)
Δρ_full[find_unstable]
Θₗ_full[find_unstable]

##

Θₗ_hist = fit(Histogram, Θₗ_full, Θ_lower_range)

plot(Θₗ_hist)

Θₗ_binidx = StatsBase.binindex.(Ref(Θₗ_hist), Θₗ_full)

for i ∈ eachindex(Θ_lower_range)

    match_Θₗ = findall(Θₗ_binidx .== i)
    find_unstable = findall(Δρ_full[match_Θₗ] .> Δρ_thres[i])
    println("$(i): $(length(find_unstable))")

end
## Take index = 21 => Θ_lower_range[21] = 0.543939393939394
lower_idx = 25
match_bins = findall(Θₗ_binidx .== lower_idx)
find_cabbeling_unstable = findall(Δρ_full[match_bins] .> Δρ_thres[lower_idx])

profile_idx = 7
Δρ_full[match_bins][find_cabbeling_unstable[profile_idx]]

long = lons_full[match_bins][find_cabbeling_unstable[profile_idx]]
lat = lats_full[match_bins][find_cabbeling_unstable[profile_idx]]
day = dates_full[match_bins][find_cabbeling_unstable[profile_idx]]
Θₗ = Θₗ_full[match_bins][find_cabbeling_unstable[profile_idx]]
Θᵤ = Θᵤ_full[match_bins][find_cabbeling_unstable[profile_idx]]
