## Make the main schematic for this work
using .VerticalProflieStability
using ColorSchemes

##
density_grad = get(ColorSchemes.dense, range(0.25, 1, length = 3))
fig = Figure(resolution = (1200, 700))
##
Θ = range(-2, 3, 1000)
Θ_grid = Θ' .* ones(length(Θ))
S = range(34.4, 34.925, 1000)
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
S_linear = range(34.514, S[end]; length = length(iso_S))
## stability schematic plot, reverse order beceause of variables
ax2 = Axis(fig[1, 2],
        title = "(b) Stability schematic",
        xlabel = "Absolute salinity (g/kg)",
        xticksvisible = false,
        xticklabelsvisible = false,
        xgridvisible = false,
        #ylabel = "Conservative temperature (∘C)",
        yticksvisible = false,
        yticklabelsvisible = false,
        ygridvisible = false
        )

# contour!(ax2, S, Θ, ρ;
#         label = L"Isopycnal through $(S^{*},~\Theta^{*})$",
#         color = density_grad[2],
#         linewidth = 2,
#         levels = [lower_isopycnal])
lines!(ax2, iso_S, iso_Θ; color = density_grad[2], linewidth = 2,
        label = L"Isopycnal through $(S^{*},~\Theta^{*})$")

# for the bands below
iso_lin_Θ = @. Θₗ + m * (iso_S - Sₗ)
Θ_linear = @. Θₗ + m * (S_linear - Sₗ)
lines!(ax2, S_linear, Θ_linear; color = density_grad[1], linewidth = 2,
       label = L"Linearised density at $(S^{*},~\Theta^{*})$")
scatter!(ax2, [Sₗ], [Θₗ]; color = density_grad[2])

# bands
## cabbeling, between the curves
band!(ax2, iso_S, iso_Θ, iso_lin_Θ; color = (density_grad[2], 0.25))
# linear extreme to start of curves
fill_S = range(S_linear[1], iso_S[1]; length = 10)
fill_Θ = fill(minimum(Θ), length(fill_S))
fill_Θ_linear = @. Θₗ + m * (fill_S - Sₗ)
band!(ax2, fill_S, fill_Θ, fill_Θ_linear; color = (density_grad[2], 0.25))
## stable
S_stable = range(minimum(S), S_linear[1]; length = length(iso_S))
Θ_stable_upper = fill(maximum(Θ_linear), length(iso_S))
Θ_stable_lower = fill(minimum(Θ), length(iso_S))
band!(ax2, S_stable, Θ_stable_lower, Θ_stable_upper; color = (density_grad[1], 0.25))
S_stable_2 = range(S_linear[1], maximum(S); length = length(iso_S))
band!(ax2, S_linear, Θ_linear, Θ_stable_upper; color = (density_grad[1], 0.25))
## unstable
Θ_unstable_fill = fill(Θ[1], length(iso_S))
band!(ax2, iso_S, Θ_unstable_fill, iso_Θ; color = (density_grad[end], 0.25))

## range bars, not sure if these add anything
# ΔΘ_stem = 2
# rangebars!(ax2, [Sₗ], [Θₗ - ΔΘ_stem], [Θₗ + ΔΘ_stem]; whiskerwidth = 10, color = :orange)
# fig

# points on linear density displaced by ± ΔΘ_stem
ΔΘ_stem = 2
S_ΔΘ = [Sₗ - (αₗ / βₗ) * ΔΘ_stem, Sₗ + (αₗ / βₗ) * ΔΘ_stem]
ΔΘ_vec = [Θₗ - ΔΘ_stem, Θₗ + ΔΘ_stem]
scatter!(ax2, S_ΔΘ, ΔΘ_vec; color = density_grad[1])

# text
## label blue dot as star
water_parcel_pos = (Sₗ+0.005, Θₗ)
wm_star = L"\left(S^{*},~\Theta^{*}\right)"
text!(ax2, water_parcel_pos, text = wm_star, fontsize = 20, align = (:left, :center),
     color = density_grad[2])
## + ΔΘ
upper_ΔΘ = L"\left(S^{ΔΘ},\Theta^{*}+\Delta\Theta\right)"
text!(ax2, S_ΔΘ[2]-0.005, ΔΘ_vec[2], text = upper_ΔΘ,
      fontsize = 20, align = (:right, :center), color = density_grad[1])
lower_ΔΘ = L"\left(S^{-\Delta\Theta},\Theta^{*}-\Delta\Theta\right)"
text!(ax2, S_ΔΘ[1]-0.005, ΔΘ_vec[1], text = lower_ΔΘ,
      fontsize = 20, align = (:right, :center), color = density_grad[1])
## Stability
text!(ax2, S[150], Θ[750], text = "Statically stable,\nstable to cabbeling",
      color = density_grad[1])
text!(ax2, S[900], Θ[200], text = "Statically unstable", color = density_grad[end],
      align = (:right, :center))
text!(ax2, S[80], Θ[235], text = "Statically stable,\nunstable to cabbeling",
      color = density_grad[2])
arrows!(ax2, [S[200]], [Θ[220]], [S[400]-S[200]], [Θ[100] - Θ[200]];
        lengthscale = 0.77, arrowcolor = density_grad[2], linecolor = density_grad[2])
arrows!(ax2, [S[200]], [Θ[320]], [S[800]-S[200]], [Θ[700] - Θ[200]];
        lengthscale = 1.07, arrowcolor = density_grad[2], linecolor = density_grad[2])
## midpoint label for ΔΘ, if using rangebars
# text!(ax2, water_parcel_pos[1], ΔΘ_vec[1] / 2, text = L"\Delta\Theta")
# text!(ax2, Sₗ-0.005, ΔΘ_vec[2] / 2, text = L"\Delta\Theta", align = (:right, :center))
#fig
axislegend(ax2, position = :lt)
#linkyaxes!(ax2, ax1)

## Cabbeling from mixing
ax1 = Axis(fig[1, 1],
        title = "(a) Mixing and the non-linear EOS",
        xlabel = "Absolute salinity (g/kg)",
        xticksvisible = false,
        xticklabelsvisible = false,
        xgridvisible = false,
        ylabel = "Conservative temperature (∘C)",
        yticksvisible = false,
        yticklabelsvisible = false,
        ygridvisible = false
        )
right_shift = 250
# contour!(ax1, S[right_shift:end], Θ, ρ[right_shift:end, :];
#         color = density_grad[1],
#         linewidth = 2,
#         levels = [lower_isopycnal],
#         label = "Isopycnal")
lines!(ax1, iso_S, iso_Θ; color = density_grad[1], linewidth = 2, label = "Isopycnal")
Θᵤ = 2.5
find_Θᵤ = findfirst(Θ .≥ Θᵤ)
find_Sᵤ = findfirst(ρ[:, find_Θᵤ] .≥ lower_isopycnal)
Sᵤ = S[find_Sᵤ]
scatter!(ax1, [Sᵤ], [Θᵤ]; color = :red, label = "Warm/salty water parcel")
Θₗ_1 = -1.5
find_Θₗ = findfirst(Θ .≥ Θₗ_1)
find_Sₗ = findfirst(ρ[:, find_Θₗ] .≥ lower_isopycnal)
Sₗ_1 = S[find_Sₗ]
scatter!(ax1, [Sₗ_1], [Θₗ_1]; color = :blue, label = "Cold/fresh water parcel")
lines!(ax1, [Sₗ_1, Sᵤ], [Θₗ_1, Θᵤ]; color = density_grad[end],
      linestyle = :dash)
#arrows!(ax1, [Sₗ], [Θₗ], [m * (Sₗ + 0.1)], [(1/m) * (Θₗ - 1)]; lengthscale = 0.0001)
# shading
## stable
S_stable = range(minimum(S[right_shift:end]), iso_S[1]; length = length(iso_S))
Θ_stable_upper = fill(maximum(Θ), length(iso_S))
Θ_stable_lower = fill(minimum(Θ), length(iso_S))
band!(ax1, S_stable, Θ_stable_lower, Θ_stable_upper; color = (density_grad[1], 0.25))
S_stable_2 = range(S_linear[1], maximum(S); length = length(iso_S))
band!(ax1, iso_S, iso_Θ, Θ_stable_upper; color = (density_grad[1], 0.25))
## unstable
Θ_unstable_fill = fill(Θ[1], length(iso_S))
band!(ax1, iso_S, Θ_unstable_fill, iso_Θ; color = (density_grad[end], 0.25))
#text
text!(ax1, S[350], Θ[710], text = "Statically stable",
      color = density_grad[1])
text!(ax1, S[900], Θ[175], text = "Statically unstable", color = density_grad[end],
      align = (:right, :center))
axislegend(ax1, position = :lt)

##
colsize!(fig.layout, 1, Auto(0.5))
fig
#save(joinpath(PLOTDIR, "schematic.png"), fig)
