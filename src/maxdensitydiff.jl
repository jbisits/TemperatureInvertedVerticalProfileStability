# Module with functions to calculate the maximum density difference between two levels of a
# vertical profile.
module MaximumDensityDifference

using Rasters, OceanRasterConversions, GibbsSeaWater

export series_max_Δρ, series2vec, get_lats

"""
    function Δρ_cab(Sₐ::Vector, Θ::Vector, p::Vector)
Cabbeling density difference is between two vertical levels of a profile. Defined as
Δρᶜ = ρₗ(βₗ(Sᵤ - Sₗ) - αₗ(Θᵤ - Θₗ)) where ₗ is the lower level and ᵤ is the upper level.
"""
function Δρ_cab(Sₐ::Vector, Θ::Vector, p::Vector)

    ρₗ = gsw_rho(Sₐ[2], Θ[2], p[2])
    αₗ, βₗ = gsw_alpha(Sₐ[2], Θ[2], p[2]), gsw_beta(Sₐ[2], Θ[2], p[2])
    ΔSₐ = Sₐ[1] - Sₐ[2]
    ΔΘ_ = Θ[1] - Θ[2]

    return ρₗ * (βₗ * ΔSₐ - αₗ * ΔΘ_)

end

"""
    function Δρ_static(Sₐ::Vector, Θ::Vector, p::Vector)
Static density difference between two vertical levels of a profile. Defined as Δρ = ρᵤ - ρₗ
where where ₗ is the lower level and ᵤ is the upper level using midpoint pressure.
"""
function Δρ_static(Sₐ::Vector, Θ::Vector, p::Vector)

    p̄ = 0.5 * (p[1] + p[2])
    ρᵤ, ρₗ = gsw_rho(Sₐ[1], Θ[1], p̄), gsw_rho(Sₐ[2], Θ[2], p̄)

    return ρᵤ - ρₗ

end

"""
    function profiles_Δρ_max(Sₐ::Vector, Θ::Vector, p::Vector, ΔΘ_thres::Float64))
Compute the maximum cabbeling density difference between all pairs of vertical levels from
a profile. The density differences are computed using `Δρ_cab` and `Δρ_static`.
"""
function Δρ_max(Sₐ::Vector, Θ::Vector, p::Vector, ΔΘ_thres::Float64)

    first_above_ΔΘ_thres = first_ΔΘ_thres(Θ, ΔΘ_thres)
    Δρᶜ_max, Δρˢ_max, upper_level, lower_level = missing, missing, missing, missing
    if !ismissing(first_above_ΔΘ_thres)

        upper_level, lower_level = first_above_ΔΘ_thres
        Δρᶜ_max = Δρ_cab([Sₐ[upper_level], Sₐ[lower_level]],
                         [Θ[upper_level], Θ[lower_level]],
                         [p[upper_level], p[lower_level]])
        Δρˢ_max = Δρ_static([Sₐ[upper_level], Sₐ[lower_level]],
                            [Θ[upper_level], Θ[lower_level]],
                            [p[upper_level], p[lower_level]])

        for i ∈ upper_level:length(Sₐ[1:end-2])

            for j ∈ (i + 1):length(Sₐ)

                if check_ΔΘ_thres([Θ[i], Θ[j]], ΔΘ_thres)

                    Δρˢ_test = Δρ_static([Sₐ[i], Sₐ[j]], [Θ[i], Θ[j]], [p[i], p[j]])
                    if Δρˢ_test > Δρˢ_max
                        Δρᶜ_max = Δρ_cab([Sₐ[i], Sₐ[j]], [Θ[i], Θ[j]], [p[i], p[j]])
                        Δρˢ_max = Δρˢ_test
                        upper_level = i
                        lower_level = j
                    end

                end

            end

        end

    end

    return (Δρᶜ_max = Δρᶜ_max,
            Δρˢ_max = Δρˢ_max,
            upper_level = upper_level,
            lower_level = lower_level)

end

"""
    function profiles_Δρ_max(stack::RasterStack, zdepth::Float64, ΔΘ_thres::Float64)
Compute the static and cabbeling density difference from data that is saved as a
`RasterStack`. This is used to compute the density difference for data from ECCO.
"""
function profiles_Δρ_max(stack::RasterStack, zdepth::Float64, ΔΘ_thres::Float64)

    # Get the `Raster`s and dimensions
    Sₐ, Θ, p = stack[:Sₐ], stack[:Θ], stack[:p]
    lons, lats, z = dims(Sₐ, X), dims(Sₐ, Y), dims(Sₐ, Z)
    z_find = findall(z .>= zdepth)

    Δρ_cab_max = similar(Array(Θ[Z(1)]))
    Δρ_static_max = similar(Array(Θ[Z(1)]))
    upper_level_idx = similar(Array(Θ[Z(1)]), Union{Missing, Int64})
    lower_level_idx = similar(Array(Θ[Z(1)]), Union{Missing, Int64})

    for i ∈ eachindex(lons), j ∈ eachindex(lats)

        #Sₐ at `lon`, `lat` to `zdepth`
        Sₐ_ll = collect(skipmissing(Sₐ[i, j, z_find, :]))
        Θ_ll = collect(skipmissing(Θ[i, j, z_find, :]))
        p_ll = collect(skipmissing(p[i, j, z_find, :]))

        if isempty(Sₐ_ll)
            Δρ_cab_max[i, j, :] .= missing
            Δρ_static_max[i, j, :] .= missing
            upper_level_idx[i, j, :] .= missing
            lower_level_idx[i, j, :] .= missing
        else
            Δρ_max_res = Δρ_max(Sₐ_ll, Θ_ll, p_ll, ΔΘ_thres)
            Δρ_cab_max[i, j, :] .= Δρ_max_res.Δρᶜ_max
            Δρ_static_max[i, j, :] .= Δρ_max_res.Δρˢ_max
            upper_level_idx[i, j, :] .= Δρ_max_res.upper_level
            lower_level_idx[i, j, :] .= Δρ_max_res.lower_level
        end

    end

    return (Δρ_cab_max = Δρ_cab_max,
            Δρ_static_max = Δρ_static_max,
            upper_level_idx = upper_level_idx,
            lower_level_idx = lower_level_idx)

end

"""
    function series_max_Δρ(raster_series::RasterSeries; zdepth = -1000.0)
Compute the static and cabbeling density difference from data that is saved as a
`RasterSeries`.
"""
function series_max_Δρ(raster_series::RasterSeries, ΔΘ_thres::Float64; zdepth = -1000.0)

    var_names = [:Δρ_cab, :Δρ_static, :Θᵤ, :Θₗ, :ΔΘ]
    var_mats = Array{Array}(undef, length(var_names))
    x, y, time = dims(raster_series[Ti(1)], X), dims(raster_series[Ti(1)], Y),
                 dims(raster_series[Ti(1)], Ti)
    timestamps = dims(raster_series, Ti)
    dd_rs_stacks = Vector{RasterStack}(undef, length(raster_series))

    for (i, stack) ∈ enumerate(raster_series)

        @info "Date $(timestamps[i])"

        converted_stack = convert_ocean_vars(stack, (Sₚ = :SALT, θ = :THETA))
        profile_max_res = profiles_Δρ_max(converted_stack, zdepth, ΔΘ_thres)
        var_mats[1] = profile_max_res.Δρ_cab_max
        var_mats[2] = profile_max_res.Δρ_static_max
        var_mats[3] = Θ_upper(converted_stack[:Θ], profile_max_res.upper_level_idx)
        var_mats[4] = Θ_lower(converted_stack[:Θ], profile_max_res.lower_level_idx)
        var_mats[5] = ΔΘ(converted_stack[:Θ], profile_max_res.upper_level_idx,
                                     profile_max_res.lower_level_idx)
        rs = [Raster(var_mats[j], (x, y, time); name = var_names[j])
                for j ∈ eachindex(var_mats)]
        dd_rs_stacks[i] = RasterStack(rs...)

    end

    return RasterSeries(dd_rs_stacks, timestamps)

end

"""
    function first_ΔΘ_thres(Θ::Vector, ΔΘ_thres::Float64)
Find the first pair of vertically separated levels of a water column that have an absolute
temperature difference greater than `ΔΘ_thres`.
"""
function first_ΔΘ_thres(Θ::Vector, ΔΘ_thres::Float64)

    upper_level = 1
    find_ΔΘ = findfirst(abs.(Θ[upper_level] .- Θ) .> ΔΘ_thres)

    upper_level += 1
    while isnothing(find_ΔΘ) && upper_level < length(Θ)

        find_ΔΘ = findfirst(abs.(Θ[upper_level] .- Θ) .> ΔΘ_thres)
        upper_level += 1

    end

    upper_level -= 1
    first_ΔΘ = isnothing(find_ΔΘ) ? missing : [upper_level, find_ΔΘ]

    return first_ΔΘ

end

"""
    function check_ΔΘ_thres(Θ::Vector, ΔΘ_thres::Float64)
Check that two vertically spaced levels have an absolute temperature difference greater than
`ΔΘ_thres`.
"""
function check_ΔΘ_thres(Θ::Vector, ΔΘ_thres::Float64)

    above_thres = false

    if abs(Θ[1] - Θ[2]) > ΔΘ_thres

        above_thres = true

    end

    return above_thres

end

"""
    function Θ_lower(Θ::Raster, lower_idx::Array)
Find temperature of lower water mass used for the maximum density difference calculation.
"""
function Θ_lower(Θ::Raster, lower_idx::Array)

    lons, lats = dims(Θ, X), dims(Θ, Y)
    Θₗ_mat = similar(Array(Θ[Z(1)]))
    for i ∈ eachindex(lons), j ∈ eachindex(lats)

        idx = lower_idx[i, j, 1]
        if !ismissing(idx)
            Θₗ_mat[i, j, 1] = Θ[i, j, idx, 1]
        end

    end

    return Θₗ_mat
end

"""
    function Θ_upper(Θ::Raster, upper_idx::Array)
Find temperature of upper water mass used for the maximum density difference calculation.
"""
function Θ_upper(Θ::Raster, upper_idx::Array)

    lons, lats = dims(Θ, X), dims(Θ, Y)
    Θᵤ_mat = similar(Array(Θ[Z(1)]))
    for i ∈ eachindex(lons), j ∈ eachindex(lats)

        idx = upper_idx[i, j, 1]
        if !ismissing(idx)
            Θᵤ_mat[i, j, 1] = Θ[i, j, idx, 1]
        end

    end

    return Θᵤ_mat

end

"""
    function ΔΘ(Θ::Raster, upper_idx::Array, lower_idx::Array)
Compute the temperature difference between the two levels at which the `max_density_angle`
is computed.
"""
function ΔΘ(Θ::Raster, upper_idx::Array, lower_idx::Array)

    ΔΘ_mat = Θ_upper(Θ, upper_idx) - Θ_lower(Θ, lower_idx)

    return ΔΘ_mat

end

"""
    function Δz(z::Vector, upper_idx::Array, lower_idx::Array)
Compute the depth difference between the two levels at which the `max_density_angle` is
calculated.
"""
function Δz(Θ::Raster, upper_idx::Array, lower_idx::Array)

    lons, lats, z = dims(Θ, X), dims(Θ, Y), dims(Θ, Z)
    Δz_mat = similar(Array(Θ[Z(1), Ti(1)]))
    for i ∈ eachindex(lons), j ∈ eachindex(lats)

        if !ismissing(upper_idx[i, j])
            Δz_mat[i, j] = abs(z[upper_idx[i, j]] - z[lower_idx[i, j]])
        end

    end

    return Δz_mat

end

"""
    function series2vec(series::RasterSeries, var::Symbol)
Convert a variable `var` from a `RasterSeries` into a `Vector` for plotting.
"""
function series2vec(series::RasterSeries, var::Symbol)

    var_vals_vec = []
    for i ∈ eachindex(series)
        var_vals = reshape(series[i][var], :)[:]
        find = @. !ismissing(var_vals)
        push!(var_vals_vec, var_vals[find])
    end

    return vcat(var_vals_vec...)

end

"""
    function get_lats(series::RasterSeries, var::Symbol)
Get the latitude at the non-missing values of a variabla `var` from a `RasterSeries`
"""
function get_lats(series::RasterSeries, var::Symbol)

    lon, lat = lookup(series[1][var], X), lookup(series[1][var], Y)
    lats_vec = repeat(lat, outer = length(lon))
    lats_vec_nm = []
    for i ∈ eachindex(series)
        var_vals = reshape(series[i][var], :)[:]
        find = @. !ismissing(var_vals)
        push!(lats_vec_nm, lats_vec[find])
    end

    return vcat(lats_vec_nm...)

end

end # module
