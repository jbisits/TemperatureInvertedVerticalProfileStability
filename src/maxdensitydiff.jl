# Module with functions to calculate the maximum density difference between two levels of a
# vertical profile.
module MaximumDensityDifference

using Rasters, OceanRasterConversions, GibbsSeaWater, MAT, NCDatasets, JLD2

export series_max_Δρ, series2vec, get_lats, argo_max_Δρ, goship_max_Δρ, en4_max_Δρ,
       group_ecco_ΔΘ, group_goship_ΔΘ

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
    function Δρ_max(Sₐ::Vector, Θ::Vector, p::Vector,
                    ΔΘ_thres::Union{Float64, Vector{Float64}})
Compute the maximum cabbeling density difference between all pairs of vertical levels from
a profile. The density differences are computed using `Δρ_cab` and `Δρ_static`.
"""
function Δρ_max(Sₐ::Vector, Θ::Vector, p::Vector, ΔΘ_thres::Union{Float64, Vector{Float64}})

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
    function profiles_Δρ_max(stack::RasterStack, zdepth::Float64
                             ΔΘ_thres::Union{Float64, Vector{Float64}})
Compute the static and cabbeling density difference from data that is saved as a
`RasterStack`. This is used to compute the density difference for data from ECCO.
"""
function profiles_Δρ_max(stack::RasterStack, zdepth::Float64,
                         ΔΘ_thres::Union{Float64, Vector{Float64}})

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
    function series_max_Δρ(raster_series::RasterSeries,
                           ΔΘ_thres::Union{Float64, Vector{Float64}};
                           zdepth = -1000.0)
    function series_max_Δρ(raster_series::RasterSeries,
                           ΔΘ_thres::Union{Float64, Vector{Float64}},
                           savepath::AbstractString; zdepth = -1000.0, filetype = ".nc")
Compute the static and cabbeling density difference from data that is saved as a
`RasterSeries`. Adding the argument `savepath` will save each `RasterStack` from the input
`raster_series` as a separate file. By default filetype is `.nc` though `.jld2` should also
work, just require a different workflow to open and look at etc.
"""
function series_max_Δρ(raster_series::RasterSeries,
                       ΔΘ_thres::Union{Float64, Vector{Float64}}; zdepth = -1000.0)

    var_names = [:Δρ_cab, :Δρ_static, :Θᵤ, :Θₗ, :ΔΘ, :Sₗ, :pᵤ, :pₗ, :Δp]
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
        var_mats[6] = S_lower(converted_stack[:Sₐ], profile_max_res.lower_level_idx)
        var_mats[7] = p_upper(converted_stack[:p], profile_max_res.upper_level_idx)
        var_mats[8] = p_lower(converted_stack[:p], profile_max_res.lower_level_idx)
        var_mats[9] = Δp(converted_stack[:p], profile_max_res.upper_level_idx,
                                              profile_max_res.lower_level_idx)
        rs = [Raster(var_mats[j], (x, y, time); name = var_names[j])
                for j ∈ eachindex(var_mats)]
        dd_rs_stacks[i] = RasterStack(rs...)

    end

    return RasterSeries(dd_rs_stacks, timestamps)

end
function series_max_Δρ(raster_series::RasterSeries,
                       ΔΘ_thres::Union{Float64, Vector{Float64}},
                       savepath::AbstractString; zdepth = -1000.0, filetype = ".nc")

    var_names = [:Δρ_cab, :Δρ_static, :Θᵤ, :Θₗ, :ΔΘ, :Sₗ, :pᵤ, :pₗ, :Δp]
    var_mats = Array{Array}(undef, length(var_names))
    x, y, time = dims(raster_series[Ti(1)], X), dims(raster_series[Ti(1)], Y),
                 dims(raster_series[Ti(1)], Ti)
    timestamps = dims(raster_series, Ti)

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
        var_mats[6] = S_lower(converted_stack[:Sₐ], profile_max_res.lower_level_idx)
        var_mats[7] = p_upper(converted_stack[:p], profile_max_res.upper_level_idx)
        var_mats[8] = p_lower(converted_stack[:p], profile_max_res.lower_level_idx)
        var_mats[9] = Δp(converted_stack[:p], profile_max_res.upper_level_idx,
                                              profile_max_res.lower_level_idx)
        rs = [Raster(var_mats[j], (x, y, time); name = var_names[j])
                for j ∈ eachindex(var_mats)]
        save_stack = RasterStack(rs...)
        write(joinpath(savepath, "$(timestamps[i])_$(ΔΘ_thres)"*filetype), save_stack)

    end

    return nothing

end

"""
    function first_ΔΘ_thres(Θ::Vector, ΔΘ_thres::Float64)
    function first_ΔΘ_thres(Θ::Vector, ΔΘ_thres::Vector{Float64})
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
function first_ΔΘ_thres(Θ::Vector, ΔΘ_thres::Vector{Float64})

    upper_level = 1
    find_ΔΘ = findfirst(ΔΘ_thres[1] .< abs.(Θ[upper_level] .- Θ) .< ΔΘ_thres[2])

    upper_level += 1
    while isnothing(find_ΔΘ) && upper_level < length(Θ)

        find_ΔΘ = findfirst(ΔΘ_thres[1] .< abs.(Θ[upper_level] .- Θ) .< ΔΘ_thres[2])
        upper_level += 1

    end

    upper_level -= 1
    first_ΔΘ = isnothing(find_ΔΘ) ? missing : [upper_level, find_ΔΘ]

    return first_ΔΘ

end
"""
    function check_ΔΘ_thres(Θ::Vector, ΔΘ_thres::Float64)
    function check_ΔΘ_thres(Θ::Vector, ΔΘ_thres::Vector{Float64})
Check that two vertically spaced levels have an absolute temperature difference greater than
`ΔΘ_thres`.
"""
function check_ΔΘ_thres(Θ::Vector, ΔΘ_thres::Float64)

    above_thres = false

    if abs(Θ[1] - Θ[2]) ≥ ΔΘ_thres

        above_thres = true

    end

    return above_thres

end
function check_ΔΘ_thres(Θ::Vector, ΔΘ_thres::Vector{Float64})

    above_thres = false

    if ΔΘ_thres[1] ≤ abs(Θ[1] - Θ[2]) < ΔΘ_thres[2]

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
Compute the temperature difference between the two levels at which the maximum static density
is computed.
"""
function ΔΘ(Θ::Raster, upper_idx::Array, lower_idx::Array)

    ΔΘ_mat = Θ_upper(Θ, upper_idx) - Θ_lower(Θ, lower_idx)

    return ΔΘ_mat

end

"""
    function S_lower(Sₐ::Raster, lower_idx::Array)
Find salinity of lower water mass used for the maximum density difference calculation.
"""
function S_lower(Sₐ::Raster, lower_idx::Array)

    lons, lats = dims(Sₐ, X), dims(Sₐ, Y)
    Sₗ_mat = similar(Array(Sₐ[Z(1)]))
    for i ∈ eachindex(lons), j ∈ eachindex(lats)

        idx = lower_idx[i, j, 1]
        if !ismissing(idx)
            Sₗ_mat[i, j, 1] = Sₐ[i, j, idx, 1]
        end

    end

    return Sₗ_mat
end

"""
    function p_lower(p::Raster, lower_idx::Array)
Find pressure of lower water mass used for the maximum density difference calculation.
"""
function p_lower(p::Raster, lower_idx::Array)

    lons, lats = dims(p, X), dims(p, Y)
    pₗ_mat = similar(Array(p[Z(1)]))
    for i ∈ eachindex(lons), j ∈ eachindex(lats)

        idx = lower_idx[i, j, 1]
        if !ismissing(idx)
            pₗ_mat[i, j, 1] = p[i, j, idx, 1]
        end

    end

    return pₗ_mat
end

"""
    function p_upper(p::Raster, upper_idx::Array)
Find pressure of upper water mass used for the maximum density difference calculation.
"""
function p_upper(p::Raster, upper_idx::Array)

    lons, lats = dims(p, X), dims(p, Y)
    pᵤ_mat = similar(Array(p[Z(1)]))
    for i ∈ eachindex(lons), j ∈ eachindex(lats)

        idx = upper_idx[i, j, 1]
        if !ismissing(idx)
            pᵤ_mat[i, j, 1] = p[i, j, idx, 1]
        end

    end

    return pᵤ_mat

end
"""
    function Δp(p::Raster, upper_idx::Array, lower_idx::Array)
Compute the pressure difference between the two levels at which the maximum static density
is computed.
"""
function Δp(p::Raster, upper_idx::Array, lower_idx::Array)

    Δp_mat = abs.(p_upper(p, upper_idx) - p_lower(p, lower_idx))

    return Δp_mat

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
Get the latitude at the non-missing values of a variable `var` from a `RasterSeries`
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
function get_lons(series::RasterSeries, var::Symbol)

    lon, lat = lookup(series[1][var], X), lookup(series[1][var], Y)
    lons_vec = repeat(lon, outer = length(lat))
    lons_vec_nm = []
    for i ∈ eachindex(series)
        var_vals = reshape(series[i][var], :)[:]
        find = @. !ismissing(var_vals)
        push!(lons_vec_nm, lons_vec[find])
    end

    return vcat(lons_vec_nm...)

end

"""
    function argo_max_Δρ(data_file::AbstractString,
                         ΔΘ_thres::Union{Float64, Vector{Float64}}; max_pressure = 1000)
Calculate the maximum static density difference from the `Argo_JJASO.mat` data. The function
takes in the data file, extracts the variables, converts to TEOS-10 standard then uses
`Δρ_max` to find the maximum static and cabbeling density difference, and the indexes of the
levels used to find these density differences. As with the ECCO data we only look to a depth
(in this case the vertical coordinate is pressure) of 1000m so max_pressure is 1000dbar.
"""
function argo_max_Δρ(data_file::AbstractString, ΔΘ_thres::Union{Float64, Vector{Float64}};
                     max_pressure = 1000)

    vars = matread(data_file)
    lat = vec(vars["lat"])
    lon = vec(vars["lon"])
    Sₚ = vars["sal"]
    θ = vars["temp"]
    p = vars["press"]
    Δρˢ, Δρᶜ = similar(lat, Union{Float64, Missing}), similar(lat, Union{Float64, Missing})
    Θₗ, Θᵤ = similar(Δρˢ), similar(Δρˢ)
    pₗ, pᵤ = similar(Δρˢ), similar(Δρˢ)
    Sₗ, Sᵤ = similar(Δρˢ), similar(Δρˢ)
    for i ∈ eachindex(Δρˢ)

        if i % 1000 == 0
            @info "Profile $(i) of $(length(lat))"
        end

        p_vec = vec(p[i])
        find_1000 = findall(p_vec .≤ max_pressure)
        if !isempty(find_1000)

            p_vec = p_vec[find_1000]
            Sₚ_vec, θ_vec = vec(Sₚ[i])[find_1000], vec(θ[i])[find_1000]
            Sₐ = gsw_sa_from_sp.(Sₚ_vec, p_vec, lon[i], lat[i])
            Θ = gsw_ct_from_pt.(Sₐ, θ_vec)

            res = Δρ_max(Sₐ, Θ, p_vec, ΔΘ_thres)
            Δρᶜ[i] = res.Δρᶜ_max
            Δρˢ[i] = res.Δρˢ_max
            ul = res.upper_level
            ll = res.lower_level
            if ismissing(ul)
                Θᵤ[i], Θₗ[i] = ul, ll
                pᵤ[i], pₗ[i] = ul, ll
                Sᵤ[i], Sₗ[i] = ul, ll
            else
                Θᵤ[i], Θₗ[i] = Θ[ul], Θ[ll]
                pᵤ[i], pₗ[i] = p_vec[ul], p_vec[ll]
                Sᵤ[i], Sₗ[i] = Sₐ[ul], Sₐ[ll]
            end

        end

    end

    return Dict("Δρˢ" => Δρˢ, "Δρᶜ" => Δρᶜ, "Θᵤ" => Θᵤ, "Θₗ" => Θₗ, "pᵤ" => pᵤ, "pₗ" => pₗ,
                "Sᵤ" => Sᵤ, "Sₗ" => Sₗ, "lats" => lat)

end

"""
    function goship_max_Δρ(data_files::Vector{String}, ΔΘ_thres::Union{Float64, Vector{Float64}};
                           max_pressure = 1000,
                           var_names = ("CTDSA", "CTDCT", "CTDprs", "lonlist", "latlist"))
Calculate the maximum static density difference of data that has been collated as part of
the GO-SHIP easy ocean project. The data is required for this function is the `.mat` format.
After extracting the profile data the function `Δρ_max` is used to find the maximum static
to density difference, as well as the levels at which this is computed. The output is
returned as a `Dict`.
"""
function goship_max_Δρ(data_files::Vector{String}, ΔΘ_thres::Union{Float64, Vector{Float64}};
                        max_pressure = 1000,
                        var_names = ("CTDSA", "CTDCT", "CTDprs", "lonlist", "latlist"),
                        res_keys = ("Δρˢ", "Δρᶜ", "Θᵤ", "Θₗ", "pᵤ", "pₗ", "Sᵤ", "Sₗ",
                                    "lons", "lats"))

    res_dict = Dict{String, Vector}(key => Vector{Union{Float64, Missing}}(undef, 0)
                                    for key ∈ res_keys)

    for (k, file) ∈ enumerate(data_files)

        @info "-- File $(k) of $(length(data_files))"

        vars = matread(file)["D_reported"]
        Sₐ = typeof(vars[var_names[1]]) == Matrix{Float64} ? [vars[var_names[1]]] :
                                                              vec(vars[var_names[1]])
        Θ =  typeof(vars[var_names[2]]) == Matrix{Float64} ? [vars[var_names[2]]] :
                                                              vec(vars[var_names[2]])
        p =  typeof(vars[var_names[3]]) == Matrix{Float64} ? [vars[var_names[3]]] :
                                                              vec(vars[var_names[3]])
        lons = typeof(vars[var_names[4]]) == Vector{Float64} ? [vars[var_names[4]]] :
                                                                vec(vars[var_names[4]])
        lats = typeof(vars[var_names[5]]) == Vector{Float64} ? [vars[var_names[5]]] :
                                                                vec(vars[var_names[5]])

        for j ∈ eachindex(Sₐ)

            Δρˢ = similar(vec(p[j][1, :]), Union{Float64, Missing})
            Δρᶜ = similar(Δρˢ)
            Θₗ, Θᵤ = similar(Δρˢ), similar(Δρˢ)
            pₗ, pᵤ = similar(Δρˢ), similar(Δρˢ)
            Sₗ, Sᵤ = similar(Δρˢ), similar(Δρˢ)

            for i ∈ eachindex(Δρˢ)

                find_1000 = findall(vec(p[j][:, i]) .≤ max_pressure)

                if !isempty(find_1000)

                    Sₐ_profile = vec(Sₐ[j][find_1000, i])
                    Θ_profile  = vec(Θ[j][find_1000, i])
                    p_profile  = vec(p[j][find_1000, i])
                    res = Δρ_max(Sₐ_profile, Θ_profile, p_profile, ΔΘ_thres)
                    Δρᶜ[i] = res.Δρᶜ_max
                    Δρˢ[i] = res.Δρˢ_max
                    ul = res.upper_level
                    ll = res.lower_level
                    if ismissing(ul)
                        Θᵤ[i], Θₗ[i] = ul, ll
                        pᵤ[i], pₗ[i] = ul, ll
                        Sᵤ[i], Sₗ[i] = ul, ll
                    else
                        Sᵤ[i], Sₗ[i] = Sₐ_profile[ul], Sₐ_profile[ll]
                        Θᵤ[i], Θₗ[i] = Θ_profile[ul], Θ_profile[ll]
                        pᵤ[i], pₗ[i] = p_profile[ul], p_profile[ll]
                    end

                end

            end
            # `append!` to the results dictionary
            res = (Δρˢ, Δρᶜ, Θᵤ, Θₗ, pᵤ, pₗ, Sᵤ, Sₗ, lons[j], lats[j])
            for (d, key) ∈ enumerate(res_keys)
                append!(res_dict[key], res[d])
            end

        end

    end

    return res_dict

end

"""
    function en4_max_Δρ(data_files::Vector{String}, ΔΘ_thres::Union{Float64, Vector{Float64}};
                        max_depth = -1000)
Calculate the maximum static density difference, and other information from EN4 temperature
salinity profiles. The function takes in a `Vector` of `String`s that are data files of
`.nc` filetype. NCDatasets then reads them, aggregating over the appropriate dimension. Then
required infromation is computed and passed to `Δρ_max`. All further required information is
computed and the returned as a `Dict`.
"""
function en4_max_Δρ(data_files::Vector{String}, ΔΘ_thres::Union{Float64, Vector{Float64}};
                    max_depth = -1000)

    ds = NCDataset(data_files; aggdim = "N_PROF")
    lat = ds["LATITUDE"][:]
    lon = ds["LONGITUDE"][:]
    θ = ds["POTM_CORRECTED"][:, :]
    Sₚ = ds["PSAL_CORRECTED"][:, :]
    z = -ds["DEPH_CORRECTED"][:, :]

    Δρˢ, Δρᶜ = similar(lat, Union{Float64, Missing}), similar(lat, Union{Float64, Missing})
    Θₗ, Θᵤ = similar(Δρˢ), similar(Δρˢ)
    pₗ, pᵤ = similar(Δρˢ), similar(Δρˢ)
    Sₗ, Sᵤ = similar(Δρˢ), similar(Δρˢ)

    for i ∈ eachindex(lat)

        profile_depths = nomissing(z[:, i], NaN)
        find_1000 = findall(profile_depths .≥ max_depth)

        if !isnothing(find_1000)

            z_profile = nomissing(z[find_1000, i], NaN)
            Sₚ_profile = nomissing(Sₚ[find_1000, i], NaN)
            θ_profile = nomissing(θ[find_1000, i], NaN)
            p = gsw_p_from_z.(z_profile, lat[i])
            Sₐ = gsw_sa_from_sp.(Sₚ_profile, p, lon[i], lat[i])
            Θ = gsw_ct_from_pt.(θ_profile, Sₐ)

            if !isempty(Θ)
                res = Δρ_max(Sₐ, Θ, p, ΔΘ_thres)
                Δρᶜ[i] = res.Δρᶜ_max
                Δρˢ[i] = res.Δρˢ_max
                ul = res.upper_level
                ll = res.lower_level
                if ismissing(ul)
                    Θᵤ[i], Θₗ[i] = ul, ll
                    pᵤ[i], pₗ[i] = ul, ll
                    Sᵤ[i], Sₗ[i] = ul, ll
                else
                    Θᵤ[i], Θₗ[i] = Θ[ul], Θ[ll]
                    pᵤ[i], pₗ[i] = p[ul], p[ll]
                    Sᵤ[i], Sₗ[i] = Sₐ[ul], Sₐ[ll]
                end
            end

        end

    end

    find_nm = findall(.!ismissing.(Δρˢ))
    return Dict("Δρˢ" => Δρˢ[find_nm], "Δρᶜ" => Δρᶜ[find_nm], "Θᵤ" => Θᵤ[find_nm],
                "Θₗ" => Θₗ[find_nm], "pᵤ" => pᵤ[find_nm], "pₗ" => pₗ[find_nm],
                "Sᵤ" => Sᵤ[find_nm], "Sₗ" => Sₗ[find_nm], "lats" => lat[find_nm])
end
"""
    function group_ecco_ΔΘ(ecco_output::AbstractString)
Group the maximum static density differences from `ecco_output` by the size of the
temperature difference. The grouping is done over a ranges of temperature values.
"""
function group_ecco_ΔΘ(ecco_output::AbstractString)

    ecco = jldopen(ecco_output)
    temp_ranges = ((-1.5, -0.5), (-2.5, -1.5), (-3.5, -2.5), (-4.5, -3.5))
    ΔΘ_keys = ("ΔΘ_0.5_1.5", "ΔΘ_1.5_2.5", "ΔΘ_2.5_3.5", "ΔΘ_3.5_4.5")
    ΔΘ_dict = Dict{String, Vector}(key => Vector{Float64}(undef, 0) for key ∈ ΔΘ_keys)
    Δρ_dict = Dict{String, Vector}(key => Vector{Float64}(undef, 0) for key ∈ ΔΘ_keys)

    for k ∈ keys(ecco)

        for (i, ΔΘ) ∈ enumerate(temp_ranges)

            find_ΔΘ = findall(ΔΘ[1] .< ecco[k]["ΔΘ_vals"] .≤ ΔΘ[2])
            append!(ΔΘ_dict[ΔΘ_keys[i]], ecco[k]["ΔΘ_vals"][find_ΔΘ])
            append!(Δρ_dict[ΔΘ_keys[i]], ecco[k]["Δρˢ"][find_ΔΘ])

        end

    end

    close(ecco)

    return ΔΘ_dict, Δρ_dict

end
"""
    function group_goship_ΔΘ(goship::AbstractString)
Group the maximum static density differnce from `goship` by the size of the temperature
difference. The grouping is done over ranges of temperature values.
"""
function group_goship_ΔΘ(goship_output::AbstractString)

    goship = jldopen(goship_output)
    temp_ranges = ((-1.5, -0.5), (-2.5, -1.5), (-3.5, -2.5), (-4.5, -3.5))
    ΔΘ_keys = ("ΔΘ_0.5_1.5", "ΔΘ_1.5_2.5", "ΔΘ_2.5_3.5", "ΔΘ_3.5_4.5")
    ΔΘ_dict = Dict{String, Vector}(key => Vector{Float64}(undef, 0) for key ∈ ΔΘ_keys)
    Δρ_dict = Dict{String, Vector}(key => Vector{Float64}(undef, 0) for key ∈ ΔΘ_keys)

    for k ∈ keys(goship)

        for (i, ΔΘ) ∈ enumerate(temp_ranges)

            ΔΘ_vals = collect(skipmissing(goship[k]["Θᵤ"] .- goship[k]["Θₗ"]))
            Δρˢ = collect(skipmissing(goship[k]["Δρˢ"]))
            find_ΔΘ = findall(ΔΘ[1] .< ΔΘ_vals .≤ ΔΘ[2])
            append!(ΔΘ_dict[ΔΘ_keys[i]], ΔΘ_vals[find_ΔΘ])
            append!(Δρ_dict[ΔΘ_keys[i]], Δρˢ[find_ΔΘ])

        end

    end

    close(goship)

    return ΔΘ_dict, Δρ_dict

end

end # module
