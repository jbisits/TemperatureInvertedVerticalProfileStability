include("../src/VerticalProfileStability.jl")
using .VerticalProfileStability
using JLD2, NCDatasets, CairoMakie, Rasters, Glob, Dates
using ColorSchemes, GeoMakie, Statistics, StatsBase
using LinearAlgebra: normalize

##

const GOSHIP_DATA = joinpath(@__DIR__, "../data/analysis/goship.jld2")
const ECCO_TS_DATA_PATH = joinpath(@__DIR__, "../data/observations/ECCO_daily_mean_TS")
const ECCO_DATA_PATH = joinpath(@__DIR__, "../data/analysis/output_1.0")
const EXTRACTED_DATA_INV = joinpath(@__DIR__, "../data/analysis/ECCO_invertedΔΘ_extracted_data_fixedΔΘ_geq.jld2")
const GOSHIP_JOINED = joinpath(@__DIR__, "../data/analysis/goship_joined.jld2")
const ONEDMODEL_SIMULATIONS = ("initial_ΔΘ_0.5", "initial_ΔΘ_1.0", "initial_ΔΘ_2.0") .* "_tgrad"
const PAPER_PLOTS_PATH = joinpath(@__DIR__, "../plots/paper/")
const PAPER_PATH = joinpath(@__DIR__, "../../../Papers/PhD-paper1-CabbelingInstability")
const CURRENT_DIR = @__DIR__

# JPO publication
publication_theme = Theme(font="CMU Serif", fontsize = 20,
                          Axis=(titlesize = 22,
                                xlabelsize = 20, ylabelsize = 20,
                                xgridstyle = :dash, ygridstyle = :dash,
                                xtickalign = 0, ytickalign = 0,
                                yticksize = 7.5, xticksize = 7.5),
                          Legend=(framecolor = (:black, 0.5),
                                  backgroundcolor = (:white, 0.5),
                                  labelsize = 20),
                          Colorbar=(ticksize=16,
                                    tickalign=1,
                                    spinewidth=0.5))
set_theme!(publication_theme)

# # # Match others for thesis
# publication_theme = Theme(font="CMU Serif", fontsize = 20,
#                           Axis=(titlesize = 22,
#                                 xlabelsize = 20, ylabelsize = 20,
#                                 xgridstyle = :dash, ygridstyle = :dash,
#                                 xtickalign = 0, ytickalign = 0,
#                                 yticksize = 6.5, xticksize = 6.5),
#                           Legend=(framecolor = (:black, 0.5),
#                                   backgroundcolor = (:white, 0.5),
#                                   labelsize = 20),
#                           Colorbar=(ticksize=12,
#                                     tickalign=1,
#                                     spinewidth=0.5))
# new_theme = merge(theme_latexfonts(), publication_theme)
# set_theme!(new_theme)

############################################################################################
## Stability schematic, figure 1
############################################################################################

density_grad = get(ColorSchemes.dense, range(0.25, 1, length = 3))
haline_grad = get(ColorSchemes.haline, range(0, 1, length = 3))
#
Θ = range(-1.95, 2, 1000)
Θ_grid = Θ' .* ones(length(Θ))
S = range(34.4, 34.8, 1000)
S_grid = S .* ones(length(S))'

p_ref = 0.0 #reference pressure
ρ = gsw_rho.(S_grid, Θ_grid, p_ref)

Sₗ, Θₗ = 34.7,  0.5
lower_isopycnal = gsw_rho(Sₗ, Θₗ, p_ref)
αₗ = gsw_alpha(Sₗ, Θₗ, p_ref)
βₗ = gsw_beta(Sₗ, Θₗ, p_ref)
m = βₗ / αₗ

find_iso = findall(ρ .≈ lower_isopycnal)
iso_S = S_grid[find_iso] # salinity values for the isopycnal
iso_Θ = Θ_grid[find_iso] # Θ values for the isopycnal
S_linear = range(34.517, S[end-10]; length = length(iso_S))
# stability schematic plot, reverse order beceause of variables
fig = Figure(size = (800, 700))
ax2 = Axis(fig[1, 1],
        title = "(a) Stability schematic",
        xlabel = "Absolute salinity",
        xticksvisible = false,
        xticklabelsvisible = false,
        xgridvisible = false,
        ylabel = "Conservative temperature",
        yticksvisible = false,
        yticklabelsvisible = false,
        ygridvisible = false
        )
hidespines!(ax2)
lines!(ax2, iso_S, iso_Θ; color = density_grad[2], linewidth = 2,
        label = L"Isopycnal through $(S^{*},~\Theta^{*})$")

# for the bands below
iso_lin_Θ = @. Θₗ + m * (iso_S - Sₗ)
Θ_linear = @. Θₗ + m * (S_linear - Sₗ)
lines!(ax2, S_linear, Θ_linear; color = density_grad[1], linewidth = 2,
       label = L"Linearised density at $(S^{*},~\Theta^{*})$")

# bands
# cabbeling, between the curves
band!(ax2, iso_S, iso_Θ, iso_lin_Θ; color = (density_grad[2], 0.25))
# linear extreme to start of curves
fill_S = range(S_linear[1], iso_S[1]; length = 10)
fill_Θ = fill(minimum(Θ), length(fill_S))
fill_Θ_linear = @. Θₗ + m * (fill_S - Sₗ)
band!(ax2, fill_S, fill_Θ, fill_Θ_linear; color = (density_grad[2], 0.25))
# stable
S_stable = range(minimum(S), S_linear[1]; length = length(iso_S))
Θ_stable_upper = fill(maximum(Θ_linear), length(iso_S))
Θ_stable_lower = fill(minimum(Θ), length(iso_S))
band!(ax2, S_stable, Θ_stable_lower, Θ_stable_upper; color = (density_grad[1], 0.25))
S_stable_2 = range(S_linear[1], maximum(S); length = length(iso_S))
band!(ax2, S_linear, Θ_linear, Θ_stable_upper; color = (density_grad[1], 0.25))
# unstable
Θ_unstable_fill = fill(Θ[1], length(iso_S))
band!(ax2, iso_S, Θ_unstable_fill, iso_Θ; color = (density_grad[end], 0.25))

# Deep water mass
# scatter!(ax2, [Sₗ], [Θₗ];  color = (:red, 0.8), label = L"\text{Deep water mass}")

# points on linear density displaced by ± ΔΘ
ΔΘ = 2
S_ΔΘ = [Sₗ - (αₗ / βₗ) * ΔΘ]
ΔΘ_vec = [Θₗ - ΔΘ]
Θ_ΔΘ = findfirst(iso_Θ .>= -1.5)
cab_salinity = 0.5 * (iso_S[Θ_ΔΘ] + S_ΔΘ[1])
# bracket
bracket!(ax2, Sₗ, Θₗ, Sₗ, ΔΘ_vec[1]; text = L" ΔΘ", width = 25, rotation = 0)
hlines!(ax2, ΔΘ_vec[1], xmin = 0.315, xmax = 0.74, color = (:black, 0.5), linestyle = :dot)
scatter!(ax2, [cab_salinity] .-0.05, ΔΘ_vec; label = L"\text{Stable}", color = haline_grad[1])
scatter!(ax2, [cab_salinity], ΔΘ_vec; label = L"\text{Unstable to cabbeling}", color = haline_grad[2])
unstable_salinity = cab_salinity + 0.05
scatter!(ax2, [unstable_salinity], ΔΘ_vec; label = L"\text{Statically unstable}", color = haline_grad[3])

# text
water_parcel_pos = (Sₗ+0.005, Θₗ)
wm_star = L"\left(S^{*},~\Theta^{*}\right)"
text!(ax2, water_parcel_pos, text = wm_star, align = (:left, :center),
     color = :red)

# Mixing lines
lines!(ax2, [cab_salinity-0.05, Sₗ], [ΔΘ_vec[1], Θₗ];
      color = haline_grad[1], linestyle = :dash)
lines!(ax2, [cab_salinity, Sₗ], [ΔΘ_vec[1], Θₗ];
       color = haline_grad[2], linestyle = :dash)
lines!(ax2, [cab_salinity+0.05, Sₗ], [ΔΘ_vec[1], Θₗ];
       color = haline_grad[3], linestyle = :dash)
arrows!(ax2, [S[1]], [Θ[1]], [1], [0]; lengthscale = 0.399)
arrows!(ax2, [S[1]], [Θ[1]], [0], [1]; lengthscale = 3.75)
# Deep water mass
scatter!(ax2, [Sₗ], [Θₗ];  color = (:red, 0.8), label = L"\text{Deep water mass}")
# legend
axislegend(ax2, position = (0.060, 0.955))
fig
##
save(joinpath(CURRENT_DIR, "fig1_schematic.png"), fig, px_per_unit = 8)
############################################################################################
## Density difference threshold, figure 2
############################################################################################

S_average = 35
S_range = [S_average - 4, S_average + 4]
Θ_range = range(-2, 30; length = 1000)
pᵣ = 500
ΔΘ_vals = (0.5, 1.0, 2.0, 3.0)
colours = reverse(get(ColorSchemes.thermal, range(0, 0.8; length = 4)))

fig = Figure(size = (800, 700))
ax = Axis(fig[1, 1];
          title = "Cabbeling instability threshold",
          titlegap = 15,
          xlabel = "Deeper water temperature (ΘᵒC)",
          xaxisposition = :top,
          ylabel = L"$Δρ$ (kgm⁻³)")
xlims!(ax, -2, 30)
for (i, ΔΘ) ∈ enumerate(ΔΘ_vals)

    # mean
    α_mean = @. gsw_alpha(S_average, Θ_range, pᵣ)
    β_mean = @. gsw_beta(S_average, Θ_range, pᵣ)
    S_thres_mean = @. S_average - (α_mean / β_mean) * (ΔΘ)
    Δρ_thres_mean = @. gsw_rho(S_thres_mean, Θ_range - ΔΘ, pᵣ) -
                       gsw_rho(S_average, Θ_range, pᵣ)
    lines!(ax, Θ_range, Δρ_thres_mean; color = colours[i], label = " ΔΘ = -$(ΔΘ)°C")

    # error bands
    α_lower = @. gsw_alpha(S_range[1], Θ_range, pᵣ)
    β_lower = @. gsw_beta(S_range[1], Θ_range, pᵣ)
    S_thres_lower = @. S_range[1] - (α_lower / β_lower) * ΔΘ
    Δρ_thres_lower = @. gsw_rho(S_thres_lower, Θ_range - ΔΘ, pᵣ) -
                        gsw_rho(S_range[1], Θ_range, pᵣ)

    α_upper = @. gsw_alpha(S_range[2], Θ_range, pᵣ)
    β_upper = @. gsw_beta(S_range[2], Θ_range, pᵣ)
    S_thres_upper = @. S_range[2] - (α_upper / β_upper) * ΔΘ
    Δρ_thres_upper = @. gsw_rho(S_thres_upper, Θ_range - ΔΘ, pᵣ) -
                        gsw_rho(S_range[2], Θ_range, pᵣ)

    band!(ax, Θ_range, Δρ_thres_lower, Δρ_thres_upper; color = (colours[i], 0.2))
end
Legend(fig[2, 1], ax, "Cabbeling instability threshold for", orientation = :horizontal)
line = hlines!(ax, 0; linestyle = :dash, color = :black)
axislegend(ax, [line], ["Static instability"], position = :rb)
fig
##
save(joinpath(CURRENT_DIR, "fig2_deltarho_threshold.png"), fig, px_per_unit = 8)
#####################################################################################
## Figure 3, S-T evolution
#####################################################################################
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
##
day_idx = findfirst(day .== timestamps)
fig = Figure(size = (1500, 600))
ax = [Axis(fig[1, i], xlabel = "Practical salinity (psu)", xticklabelrotation = π / 4,
           ylabel = "Potential temperature (°C)",  width = 270, height = 420) for i ∈ 1:3]
Slimits = begin
    rs = Raster(ECCO_files[day_idx-1], name = :SALT)
    S = skipmissing(vec(rs[X(At(long)), Y(At(lat))].data))
    extrema(S)
end
for (i, d) ∈ enumerate(day_idx-1:day_idx+1)

    rs = RasterStack(ECCO_files[d], name = (:SALT, :THETA))
    S = vec(rs[:SALT][X(At(long)), Y(At(lat))].data)
    T = vec(rs[:THETA][X(At(long)), Y(At(lat))].data)
    S_plot = S[.!ismissing.(S)]
    T_plot = T[.!ismissing.(S)]
    z = lookup(rs, Z)[.!ismissing.(S)]
    z_range = findall(z .> -750)

    upper_idx = findfirst(T_plot .≥ Θᵤ)
    lower_idx = isnothing(findfirst(T_plot .≥ Θₗ)) ? findmax(T_plot)[2] : findfirst(T_plot .≥ Θₗ)

    zᵤ = z[upper_idx]
    zₗ = z[lower_idx]
    ẑ = -(zᵤ + zₗ) / 2

    ρᵤ = gsw_rho(S[upper_idx], T[upper_idx], ẑ)
    ρₗ = gsw_rho(S[lower_idx], T[lower_idx], ẑ)
    ρᵤ - ρₗ

    S_linear_range = range(Slimits..., length = 50)
    αₗ = gsw_alpha(S[lower_idx], T[lower_idx], ẑ)
    βₗ = gsw_beta(S[lower_idx], T[lower_idx], ẑ)
    Θ_linear = @. T[lower_idx] + (βₗ / αₗ) * (S_linear_range - S[lower_idx])
    lines!(ax[i], S_linear_range, Θ_linear, label = L"Tangent to isopycnal at $(S^{*}, Θ^{*})$", color = :red, linestyle = :dash)

    S_range = range(Slimits..., length = 1000)
    T_range = range(extrema(skipmissing(T))..., length = 1000)

    S_grid = S_range .* ones(length(S_range))'
    T_grid = T_range' .* ones(length(T_range))

    ρ = gsw_rho.(S_grid, T_grid, ẑ)

    contour!(ax[i], S_range, T_range, ρ, levels = [ρₗ], label = L"Isopycnal at $(S^{*}, Θ^{*})$", color = :black, linestyle = :dot)

    cmap = :viridis
    plt = scatterlines!(ax[i], S_plot[z_range], T_plot[z_range], color = z[z_range],
                       markercolormap = cmap, colormap = cmap, markersize = 7,
                       label = L"ECCO profile$$")
    text!(ax[i], 34.895, 1.15, text = "$(timestamps[d])", font = :bold, align = (:center, :center))
    # ax[i].title = "$(timestamps[d])"

    upper_idx = findfirst(T_plot .≥ Θᵤ)
    lower_idx = isnothing(findfirst(T_plot .≥ Θₗ)) ? findmax(T_plot)[2] : findfirst(T_plot .≥ Θₗ)
    if i != 3
        scatter!(ax[i], [S[upper_idx]], [T[upper_idx]], color = :blue, marker = :cross,
                label = L"Upper water$$")
    end
    scatter!(ax[i], [S[lower_idx]], [ T[lower_idx]], color = :red, marker = :cross,
            label = L"(S^{*}, Θ^{*})")

    if d > day_idx-1
        # Horizontal layout
        hideydecorations!(ax[i], grid = false)
    end
    if i != 2
        # Horizontal layout
        hidexdecorations!(ax[i], grid = false, ticks = false, ticklabels = false)
    end
    if i == 3
        # Horizontal layout
        Colorbar(fig[1, 4], limits = extrema(z[z_range]), colormap = cmap, label = "z (m)")
    end
    xlims!(ax[i], 34.87, 34.974)
    ylims!(ax[i], 0, 1.2)
end
Label(fig[0, :], "ECCO Salinity-temperature profile at $(long)°E, $(lat)°N", font = :bold)
# Horizontal layout
Legend(fig[2, :], ax[1], orientation  = :horizontal, nbanks = 2)
linkyaxes!(ax[1], ax[2])
linkyaxes!(ax[1], ax[3])
fig
resize_to_layout!(fig)
fig
##
save(joinpath(CURRENT_DIR, "fig3_cabbeling_unstable_profile_horizontal.png"), fig, px_per_unit = 8)
############################################################################################
## Figure 4, profile locations
############################################################################################

## Temperature inverted profile location map
output_path = ECCO_DATA_PATH
output_files = glob("*.nc", output_path)
ΔΘ_series = RasterSeries(output_files, Ti, child = Raster, name = :ΔΘ)

for rs ∈ ΔΘ_series
    map!(x -> ismissing(x) ? x : x ≤ -0.5 ? x : missing, rs, rs)
end

heatmap(ΔΘ_series[end])
length(ΔΘ_series)

# Change to 1, then can merge and sum over the time dimension.

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

rs_proportion = Raster(ΔΘ_counts, (X(lookup(ΔΘ_series[1], :X)), Y(lookup(ΔΘ_series[1], :Y))))

## Figure
fig = Figure(size = (500, 1000))
ECCO_plot = fig[1, 1] = GridLayout()
ax = GeoAxis(ECCO_plot[1, 1];
             title = "(a) ECCOv4r4",
             xlabel = "Longitude",
             xticklabelrotation = pi/4,
             ylabel = "Latitude",
             dest = "+proj=wintri",
             coastlines = true)
hidedecorations!(ax)
hm = heatmap!(ax, lookup(rs_proportion, :X), lookup(rs_proportion, :Y), ΔΘ_counts;
              colormap = :batlow)
Colorbar(fig[1, 2], hm, label = "Frequency")

# Scatter marker for unstable to cabbeling profile
scatter!(ax, [long], [lat], marker = :star5, color = :yellow, markersize = 18)
fig
## GOSHIP part of map

gd = jldopen(GOSHIP_DATA)

ΔΘ_vals = (0.5, 1.0, 2.0, 3.0)
oceans = ("atlantic", "indian", "pacific", "southern")
GOSHIP_plot = fig[2, 1] = GridLayout()
ax = GeoAxis(GOSHIP_plot[1, 1];
             title = "(b) GO-SHIP",
             xlabel = "Longitude",
             xticklabelrotation = pi/4,
             dest = "+proj=wintri",
             coastlines = true)
hidedecorations!(ax)
for key ∈ keys(gd["1.0"])

    Θᵤ = collect(skipmissing(gd["1.0"][key]["Θᵤ"]))
    Θₗ = collect(skipmissing(gd["1.0"][key]["Θₗ"]))
    lons = gd["1.0"][key]["lons"]
    lats = gd["1.0"][key]["lats"]
    pᵤ = gd["1.0"][key]["pᵤ"]
    pₗ = gd["1.0"][key]["pₗ"]
    find_inv = findall(Θᵤ .< Θₗ)
    lons_inv = lons[find_inv]
    lats_inv = lats[find_inv]
    scatter!(ax, lons_inv, lats_inv; color = :red, markersize = 3)
end
close(gd)

colsize!(fig.layout, 2, Aspect(1, 0.5))
colgap!(fig.layout, -50)
fig
##
save(joinpath(CURRENT_DIR, "fig4_profilelocation_withmarker.png"), fig, px_per_unit = 8)
############################################################################################
## Figure 5, model initial static density difference
############################################################################################
expt_ls = (:dot, :dash, :dashdot)
colours = reverse(get(ColorSchemes.thermal, range(0, 0.8; length = 4)))
density_grad = get(ColorSchemes.dense, range(0.25, 1, length = 3))
ΔΘ_thres_vals = (0.5, 1.0, 2.0)
fig = Figure(size = (800, 700))
xtickposition = -2:0
ax = Axis(fig[1, 1];
          title = "Initial static density difference at model interface",
          xlabel = "Initial ΔΘ (°C) at model interface",
          xticks = (xtickposition, string.([-2, -1, -0.5])),
          ylabel = L"$Δσ_{0}$ (kgm⁻³)")
hlines!(ax, 0; color = :black, linestyle = :dash)
for (j, sim) ∈ enumerate(ONEDMODEL_SIMULATIONS)
    sim_output = joinpath(SIM_DATADIR, sim)
    saved_ts = jldopen(joinpath(sim_output, "output_timeseries.jld2"))
    S_ts_, T_ts_ = saved_ts["S_ts"], saved_ts["T_ts"]
    num_ics = 2
    close(saved_ts)
    Θₗ_ = [T_ts_[80, 1, i] for i ∈ 1:num_ics]
    Sₗ_ = [S_ts_[80, 1, i] for i ∈ 1:num_ics]
    Θᵤ_ = [T_ts_[81, 1, i] for i ∈ 1:num_ics]
    Sᵤ_ = [S_ts_[81, 1, i] for i ∈ 1:num_ics]

    αₗ_ = gsw_alpha.(Sₗ_, Θₗ_, p_ref)
    βₗ_ = gsw_beta.(Sₗ_, Θₗ_, p_ref)

    Δρ_thres = @. gsw_rho(Sₗ_ - (αₗ_ / βₗ_) * ΔΘ_thres_vals[j],
                        Θₗ_ - ΔΘ_thres_vals[j], p_ref) -
                gsw_rho(Sₗ_, Θₗ_, p_ref)
    Δρ_static = @. gsw_rho(Sᵤ_, Θᵤ_, p_ref) - gsw_rho(Sₗ_, Θₗ_, p_ref)
    hlines!(ax, Δρ_thres; color = colours[j], linestyle = :dot,
            label = "initial ΔΘ = -$(ΔΘ_thres_vals[j])°C")
    sc = scatter!(ax, fill(xtickposition[4-j], 2), Δρ_static;
                  color = [:blue, :red], marker = [:circle, :rect])
    if j == 3
        axislegend(ax, ax, "Cabbeling instability threshold for"; position = (1, 0.1),
                   orientation = :horizontal, nbanks = 4)
    end
end
markers = [MarkerElement(marker = :circle, color = :blue),
           MarkerElement(marker = :rect, color = :red)]
labels = ["Stable to cabbeling", "Unstable to cabbeling"]
Legend(fig[2, 1], markers, labels, orientation = :horizontal)
fig
##
save(joinpath(CURRENT_DIR, "fig5_model_deltarho_threshold.png"), fig, px_per_unit = 8)
############################################################################################
## Figure 6, model density time series
############################################################################################
σ₀_fig = Figure(resolution = (1000, 1000))
z = -497.5:5:-2.5
z_range = 61:100
ax = [Axis(σ₀_fig[j, 1],
        xlabel = "Time (days)",
        ylabel = "Depth (metres)",
        limits = ((0, 60), (-200, -5))) for j ∈ 1:3]
ΔΘ_vals = (0.5, 1.0, 2.0)
letter_labels = ("(a)", "(b)", "(c)")
for (j, sim) ∈ enumerate(ONEDMODEL_SIMULATIONS)

    sim_output_ = joinpath(SIM_DATADIR, sim)
    saved_ts_ = jldopen(joinpath(sim_output_, "output_timeseries.jld2"))
    S_ts_, T_ts_, κ_ts_ = saved_ts_["S_ts"], saved_ts_["T_ts"], saved_ts_["κ_ts"]
    t = saved_ts_["t"]
    close(saved_ts_)
    σₒ_ts = gsw_sigma0.(S_ts_[z_range, :, :], T_ts_[z_range, :, :])

    sal_ics = round.(S_ts_[end, 1, end]; digits = 4)

    ax[j].title = "$(letter_labels[j]) Initial ΔΘ = -$(ΔΘ_vals[j])°C"
    if j != 3
        hidexdecorations!(ax[j], ticks = false)
    end
    reduced_z = z[z_range]
    upper = findfirst.(isapprox.(κ_ts_[z_range, i, 2], 1.0; atol = 1e-5)
                       for i ∈ eachindex(κ_ts_[1, :, 2]))
    upper_idx = findall(isnothing, upper)
    upper = map(x -> isnothing(x) ? NaN : x, upper)
    upper_plot = [reduced_z[upper[i]] for i ∈ findall(!isnan, upper)] .- 2.5
    lower = findlast.(isapprox.(κ_ts_[z_range, i, 2], 1.0; atol = 1e-5)
                       for i ∈ eachindex(κ_ts_[1, :, 2]))
    lower_idx = findall(isnothing, lower)
    lower = map(x -> isnothing(x) ? NaN : x, lower)
    lower_plot = [reduced_z[lower[i]] for i ∈ findall(!isnan, lower)] .+ 2.5
    p = Pattern('x'; width = 3, linecolor = :black)
    rapid_diff_region = vcat(collect(zip(t[findall(!isnan, lower)], lower_plot)),
                             reverse(collect(zip(t[findall(!isnan, upper)], upper_plot))))
    # Plotting
    poly!(ax[j], Point2f[rapid_diff_region...], color = p)
    hm = heatmap!(ax[j], t, z[z_range], σₒ_ts[:, :, 2]'; colormap = (:dense, 0.9))
    if j == length(ONEDMODEL_SIMULATIONS)
        Colorbar(σ₀_fig[4, 1], hm, label = L"$σ_{0}$ anomaly (kgm⁻³)", vertical = false,
                 flipaxis = false,
                 ticklabelrotation = π/4)
    end

end
σ₀_fig
##
save(joinpath(CURRENT_DIR, "fig6_model_rho_ts.png"), σ₀_fig, px_per_unit = 8)
############################################################################################
## Figure 6, ECCO pdf
############################################################################################
inv_data = jldopen(EXTRACTED_DATA_INV)
Δρ_limits = (-0.2, 0.01)
xlimits = (-1.88, 10)
bin_width = 0.0001
colours = reverse(get(ColorSchemes.thermal, range(0, 0.8; length = 4)))
area = begin
        grid_path = joinpath(@__DIR__, "../data/observations/ECCO_grid/GRID_GEOMETRY_ECCO_V4r4_latlon_0p50deg.nc")
        rs_grid = Raster(grid_path, name = :area)
        rs_grid[X(1)]
       end
full_fig = Figure(size = (1000, 1200))
# scatter
splot = full_fig[1:4, 1:2] = GridLayout()
ax_splot = Axis(splot[1, 1];
                xlabel = "Deeper water temperature (ΘᵒC)",
                xaxisposition = :top,
                title = "(a) Maximum static density difference, ECCOv4r4",
                titlegap = 15,
                ylabel = L"$Δρ$ (kgm⁻³)",
                limits = (xlimits, Δρ_limits))
for (i, key) ∈ enumerate(keys(inv_data))

    Θₗ, Δρˢ = inv_data[key]["Θₗ"], inv_data[key]["Δρˢ"]
    lats = inv_data[key]["lats"]
    find = findall(xlimits[1] .≤ Θₗ .≤ xlimits[2] .&& Δρ_limits[1] .≤ Δρˢ .≤ Δρ_limits[2])
    Θₗ, Δρˢ = Θₗ[find], Δρˢ[find]
    Θ_lower_range = inv_data[key]["Θ_lower_range"]
    Θ_lower_range = inv_data[key]["Θ_lower_range"]
    Δρ_thres = inv_data[key]["Δρ_thres"]
    ΔΘ_range = inv_data[key]["ΔΘ_range"]

    sc = scatter!(ax_splot, Θₗ, Δρˢ; color = colours[i], markersize = 1)
    lines!(ax_splot, Θ_lower_range, Δρ_thres; color = colours[i],
           label = "ΔΘ = -$(ΔΘ_range[1])ᵒC", linewidth = 2)

end
Legend(splot[2, 1], ax_splot, "Cabbeling instability threshold for", orientation = :horizontal, nbanks = 2)

# pdf
pdf_plots = full_fig[1:4, 3] = GridLayout()
ax_pdf = [Axis(pdf_plots[i, 1];
        xlabel = L"$Δρ$ (kgm⁻³)",
        xticklabelsize = 15, yticklabelsize = 15,
        yticksize = 7, xticksize = 7,
        titlesize = 20) for i ∈ 1:4]
less_thres = Vector{Float64}(undef, 4)
over_thres = Vector{Float64}(undef, 4)
letter_labels = ["(b)", "(c)", "(d)", "(e)"]

for (i, key) ∈ enumerate(keys(inv_data))

    Θₗ = inv_data[key]["Θₗ"]
    Δρˢ = collect(skipmissing(inv_data[key]["Δρˢ"]))
    lats = inv_data[key]["lats"]
    ΔΘ_range = inv_data[key]["ΔΘ_range"]
    Δρ_thres = inv_data[key]["Δρ_thres"]

    Δρ_thres_mean = mean(Δρ_thres)
    hist_edges = minimum(Δρˢ):bin_width:maximum(Δρˢ)
    area_weights_ = weights(Float32.([area[Y(At(lat))] for lat ∈ lats]))
    hist_fit = fit(Histogram, Δρˢ, area_weights_, hist_edges)
    hist_fit = normalize(hist_fit; mode = :pdf)
    plot!(ax_pdf[i], hist_fit; color = (colours[i], 0.5))
    vlines!(ax_pdf[i], Δρ_thres_mean; color = colours[i], linewidth = 2,
            label = L"$Δρ$ threshold for $ΔΘ$")
    vlines!(ax_pdf[i], 0; color = :black, linestyle = :dash)
    ax_pdf[i].title = letter_labels[i] * " PDF, ΔΘ = -$(ΔΘ_range[1])°C"

    find_thres = findall(hist_edges .≤ Δρ_thres_mean)

    # to average threshold
    less_thres[i] = sum(hist_fit.weights[find_thres] .* bin_width)
    # after average threshold
    over_thres[i] = 1 - less_thres[i]

end

for i ∈ 1:4
    if i != 4
        hidexdecorations!(ax_pdf[i], grid = false, ticks = false)
    end
    xlims!(ax_pdf[i], Δρ_limits)
end
close(inv_data)
# full_fig
#
save(joinpath(CURRENT_DIR, "fig7_ECCOpdfs.png"), full_fig, px_per_unit = 8)

############################################################################################
## Figure 8, GOSHIP pdfs
############################################################################################
gdj = jldopen(GOSHIP_JOINED)
Θₗ_lims = (-1.88, 10)
Θₗ_range = range(Θₗ_lims...; length = 100)
ΔΘ_thres = (0.5, 1.0, 2.0, 3.0)
ΔΘ_colours = reverse(get(ColorSchemes.thermal, range(0, 0.8; length = 4)))
Δρ_lims = (-0.2, 0.02)
bin_width = 0.001
full_fig = Figure(size = (1000, 1200))
# scatter
splot = full_fig[1:4, 1:2] = GridLayout()
ax_splot = Axis(splot[1, 1];
                xlabel = "Deeper water temperature (ΘᵒC)",
                xaxisposition = :top,
                title = "(a) Maximum static density difference, GOSHIP",
                titlegap = 15,
                ylabel = L"$Δρ$ (kgm⁻³)",
                limits = (Θₗ_lims, Δρ_lims))
data_count = Vector{Int64}(undef, length(ΔΘ_thres))
for (i, key) ∈ enumerate(keys(gdj))

    Θᵤ = collect(skipmissing(gdj[key]["Θᵤ"]))
    Θₗ = collect(skipmissing(gdj[key]["Θₗ"]))
    find_inv = Θᵤ .≤ Θₗ
    Θₗ = Θₗ[find_inv]
    Δρˢ = collect(skipmissing(gdj[key]["Δρˢ"]))[find_inv]
    data_count[i] = length(Δρˢ)

    scatter!(ax_splot, Θₗ, Δρˢ; color = ΔΘ_colours[i], markersize = 6)

    Sₗ_mean = mean(collect(skipmissing(gdj[key]["Sₗ"]))[find_inv])
    pₘ = 0.5 .* (collect(skipmissing(gdj[key]["pₗ"]))[find_inv] .+
                    collect(skipmissing(gdj[key]["pᵤ"]))[find_inv])
    pₘ_mean = mean(pₘ)
    αₗ = gsw_alpha.(Sₗ_mean, Θₗ_range, pₘ_mean)
    βₗ = gsw_beta.(Sₗ_mean, Θₗ_range, pₘ_mean)
    Δρ_thres = @. gsw_rho.(Sₗ_mean - (αₗ / βₗ) * ΔΘ_thres[i], Θₗ_range - ΔΘ_thres[i],
                            pₘ_mean) - gsw_rho.(Sₗ_mean, Θₗ_range, pₘ_mean)
    lines!(ax_splot, Θₗ_range, Δρ_thres; color = ΔΘ_colours[i], linewidth = 2,
            label = "ΔΘ = -$(ΔΘ_thres[i])ᵒC")
end
full_fig
Legend(splot[2, 1], ax_splot, "Cabbeling instability threshold for", orientation = :horizontal, nbanks = 2)
full_fig
data_count

# pdf
pdf_plots = full_fig[1:4, 3] = GridLayout()
ax_pdf = [Axis(pdf_plots[i, 1],
        titlesize = 20,
        xticklabelsize = 15, yticklabelsize = 15,
        yticksize = 7, xticksize = 7,
        xlabel = L"$Δρ$ (kgm⁻³)")
      for i ∈ 1:4]
less_thres = Vector{Float64}(undef, 4)
over_thres = Vector{Float64}(undef, 4)
data_count2 = Vector{Int64}(undef, length(ΔΘ_thres))
letter_labels = ["(b)", "(c)", "(d)", "(e)"]

for (i, key) ∈ enumerate(keys(gdj))
    Θₗ = collect(skipmissing(gdj[key]["Θₗ"]))
    Θᵤ = collect(skipmissing(gdj[key]["Θᵤ"]))
    find_inverted = Θᵤ .< Θₗ
    Δρˢ = collect(skipmissing(gdj[key]["Δρˢ"]))[find_inverted]
    data_count2[i] = length(Δρˢ)
    Sₗ_mean = mean(collect(skipmissing(gdj[key]["Sₗ"]))[find_inverted])
    pₘ = 0.5 .* (collect(skipmissing(gdj[key]["pₗ"]))[find_inverted] .+
                collect(skipmissing(gdj[key]["pᵤ"]))[find_inverted])
    pₘ_mean = mean(pₘ)
    αₗ = gsw_alpha.(Sₗ_mean, Θₗ_range, pₘ_mean)
    βₗ = gsw_beta.(Sₗ_mean, Θₗ_range, pₘ_mean)
    Δρ_thres = @. gsw_rho.(Sₗ_mean - (αₗ / βₗ) * ΔΘ_thres[i], Θₗ_range - ΔΘ_thres[i],
                        pₘ_mean) - gsw_rho.(Sₗ_mean, Θₗ_range, pₘ_mean)
    Δρ_thres_mean = mean(Δρ_thres)
    hist_edges = minimum(Δρˢ):bin_width:maximum(Δρˢ)

    hist!(ax_pdf[i], Δρˢ; bins = hist_edges, normalization = :pdf,
          color = (ΔΘ_colours[i], 0.5))
    vlines!(ax_pdf[i], Δρ_thres_mean; color = ΔΘ_colours[i], linewidth = 2,
            label = L"$Δρ$ threshold for ΔΘ")
    vlines!(ax_pdf[i], 0; color = :black, linestyle = :dash,
            label = "Static stability threshold")
    ax_pdf[i].title = letter_labels[i] * " PDF, ΔΘ = -$(ΔΘ_thres[i])°C"
    hist_fit = fit(Histogram, Δρˢ, hist_edges)
    hist_fit = normalize(hist_fit; mode = :pdf)

    find_thres = findall(hist_edges .≤ Δρ_thres_mean)

    # to average threshold
    less_thres[i] = sum(hist_fit.weights[find_thres] .* bin_width)
    # after average threshold
    over_thres[i] = 1 - less_thres[i]

end
full_fig
pdf_lims = (-0.2, 0.01)
for i ∈ 1:4
    if i != 4
        hidexdecorations!(ax_pdf[i], grid = false, ticks = false)
    end
    xlims!(ax_pdf[i], pdf_lims)
end
less_thres
over_thres
full_fig
data_count == data_count2
data_count
save(joinpath(CURRENT_DIR, "fig8_GOSHIPpdfs.png"), full_fig, px_per_unit = 8)
close(gdj)
##

############################################################################################
## Figure 9, Probability againset ΔΘ
############################################################################################
Δρ_limits = (-0.4, 0.01)
Δρ_val = -0.04
bin_width = 0.0001
colours = reverse(get(ColorSchemes.thermal, range(0, 0.8; length = 4)))
fig = Figure(size = (800, 1000))

#
inv_data = jldopen(EXTRACTED_DATA_INV)
area = begin
        grid_path = joinpath(@__DIR__, "../data/observations/ECCO_grid/GRID_GEOMETRY_ECCO_V4r4_latlon_0p50deg.nc")
        rs_grid = Raster(grid_path, name = :area)
        rs_grid[X(1)]
       end

ax_pdf = Axis(fig[1, 1];
        xlabel = L"$Δρ$ (kgm⁻³)",
        title = "(a) ECCOv4r4 pdfs")

for (i, key) ∈ enumerate(keys(inv_data))

    Θₗ = inv_data[key]["Θₗ"]
    Δρˢ = collect(skipmissing(inv_data[key]["Δρˢ"]))
    lats = inv_data[key]["lats"]
    ΔΘ_range = inv_data[key]["ΔΘ_range"]
    Δρ_thres = inv_data[key]["Δρ_thres"]

    Δρ_thres_mean = mean(Δρ_thres)
    hist_edges = minimum(Δρˢ):bin_width:maximum(Δρˢ)
    area_weights_ = weights(Float32.([area[Y(At(lat))] for lat ∈ lats]))
    hist_fit = fit(Histogram, Δρˢ, area_weights_, hist_edges)
    hist_fit = normalize(hist_fit; mode = :pdf)
    plot!(ax_pdf, hist_fit; color = (colours[i], 0.8), label = " ΔΘ = -$(ΔΘ_range[1])°C")

end
vlines!(ax_pdf, Δρ_val, color = :red, linestyle = :dash, label = L"Reference $Δρ'$")
xlims!(ax_pdf, Δρ_limits)
ylims!(ax_pdf; low = 0, high = nothing)
close(inv_data)
fig
#
ΔΘ_range = (0.5, 1.0, 2.0, 3.0)
bin_width = 0.001
ax_goship = Axis(fig[1, 2];
                xlabel = L"$Δρ$ (kgm⁻³)",
                title = "(b) GOSHIP pdfs")

gdj = jldopen(GOSHIP_JOINED)
for (i, key) ∈ enumerate(keys(gdj))

    Θₗ = collect(skipmissing(gdj[key]["Θₗ"]))
    Θᵤ = collect(skipmissing(gdj[key]["Θᵤ"]))
    find_inverted = Θᵤ .< Θₗ
    Δρˢ = collect(skipmissing(gdj[key]["Δρˢ"]))[find_inverted]
    hist_edges = minimum(Δρˢ):bin_width:maximum(Δρˢ)

    hist_fit = fit(Histogram, Δρˢ, hist_edges)
    hist_fit = normalize(hist_fit; mode = :pdf)

   plot!(ax_goship, hist_fit; color = (colours[i], 0.8), label = " ΔΘ = -$(ΔΘ_range[i])°C")

end
vlines!(ax_goship, Δρ_val, color = :red, linestyle = :dash, label = L"Δρ_{r}")
xlims!(ax_goship, Δρ_limits)
ylims!(ax_goship; low = 0, high = nothing)
close(gdj)
fig
Legend(fig[2, :],  ax_goship, orientation = :horizontal)
#
ΔΘ_vals = [0.0, -1.0, -2.0, -3.0]
ax_ecdf = Axis(fig[3, :];
           title = "(c) Temperature inversion effect on stratification",
           xlabel = L"ΔΘ~(°C~)",
           xtickformat = values -> string.([-3.0, -2.0, -1.0, -0.5]),
           ylabel = L"ℙ\left(Δρ_{\mathrm{static}}^{\mathrm{max}}~>~Δρ_{\mathrm{r}}~|~ΔΘ\right)")
scatterlines!(ax_ecdf, ΔΘ_vals, 1 .- ECCO_cdf_Δρ_val; label = "ECCOv4r4")
scatterlines!(ax_ecdf, ΔΘ_vals, 1 .- GOSHIP_cdf_Δρ_val; label = "GOSHIP")
axislegend(ax_ecdf, position = :lt)
fig
save(joinpath(CURRENT_DIR, "fig9_prob_deltatheta.png"), fig, px_per_unit = 8)
