## Make the main schematic for this work
using .VerticalProflieStability
using ColorSchemes
density_grad = get(ColorSchemes.dense, range(0.25, 1, length = 3))
##
Θ = range(-2, 3, 1000)
Θ_grid = Θ' .* ones(length(Θ))
S_grid = S .* ones(length(S))'
S = range(34.4, 34.925, 1000)

p_ref = 0.0 #reference pressure
freezing_pt = gsw_ct_freezing.(S, p_ref, 1)

ρ = gsw_rho.(S_grid, Θ_grid, p_ref)

Sₗ, Θₗ = 34.7,  0.5
lower_isopycnal = gsw_rho(Sₗ, Θₗ, p_ref)

fig = Figure(resolution = (600, 600))
ax = Axis(fig[1, 1],
        title = "Stability schematic",
        xlabel = "Absolute salinity (g/kg)",
        xticksvisible = false,
        xticklabelsvisible = false,
        xgridvisible = false,
        ylabel = "Conservative temperature (∘C)",
        yticksvisible = false,
        yticklabelsvisible = false,
        ygridvisible = false,
        #aspect = AxisAspect(1.25),
        # topspinevisible = false,
        # bottomspinevisible = false,
        # rightspinevisible = false,
        # leftspinevisible = false
        )

contour!(ax, S, Θ, ρ,
        color = density_grad[end],
        linewidth = 2,
        levels = [lower_isopycnal],
        label = "Isopycnal (1027.71)")

αₗ = gsw_alpha(Sₗ, Θₗ, p_ref)
βₗ = gsw_beta(Sₗ, Θₗ, p_ref)
m = βₗ / αₗ

# for the bands below
find_iso = findall(ρ .≈ lower_isopycnal)
iso_S = S_grid[find_iso] # salinity values for the isopycnal
iso_lin_Θ = @. Θₗ + m * (iso_S - Sₗ)
iso_Θ = Θ_grid[find_iso] # Θ values for the isopycnal
S_linear = range(34.514, S[end]; length = length(iso_S))
Θ_linear = @. Θₗ + m * (S_linear - Sₗ)
lines!(ax, S_linear, Θ_linear; color = density_grad[1])
scatter!(ax, [Sₗ], [Θₗ]; color = density_grad[3])
fig

# bands
## cabbeling, between the curves
band!(ax, iso_S, iso_Θ, iso_lin_Θ; color = (density_grad[2], 0.25))
# linear extreme to start of curves
fill_S = range(S_linear[1], iso_S[1]; length = 10)
fill_Θ = fill(minimum(Θ), length(fill_S))
fill_Θ_linear = @. Θₗ + m * (fill_S - Sₗ)
band!(ax, fill_S, fill_Θ, fill_Θ_linear; color = (density_grad[2], 0.25))
## stable
S_stable = range(minimum(S), S_linear[1]; length = length(iso_S))
Θ_stable_upper = fill(maximum(Θ_linear), length(iso_S))
Θ_stable_lower = fill(minimum(Θ), length(iso_S))
band!(ax, S_stable, Θ_stable_lower, Θ_stable_upper; color = (density_grad[1], 0.25))
S_stable_2 = range(S_linear[1], maximum(S); length = length(iso_S))
band!(ax, S_linear, Θ_linear, Θ_stable_upper; color = (density_grad[1], 0.25))
## unstable
Θ_unstable_fill = fill(Θ[1], length(S_unstable))
band!(ax, iso_S, Θ_unstable_fill, iso_Θ; color = (density_grad[end], 0.25))
fig

## range bars, not sure if these add anything
# ΔΘ_stem = 2
# rangebars!(ax, [Sₗ], [Θₗ - ΔΘ_stem], [Θₗ + ΔΘ_stem]; whiskerwidth = 10, color = :orange)
# fig

# points on linear density displaced by ± ΔΘ_stem
S_ΔΘ = [Sₗ - (αₗ / βₗ) * ΔΘ_stem, Sₗ + (αₗ / βₗ) * ΔΘ_stem]
ΔΘ_vec = [Θₗ - ΔΘ_stem, Θₗ + ΔΘ_stem]
scatter!(ax, S_ΔΘ, ΔΘ_vec; color = density_grad[1])
fig

# text
## label blue dot as star
water_parcel_pos = (Sₗ+0.005, Θₗ)
wm_star = L"\left(S^{*},~\Theta^{*}\right)"
text!(ax, water_parcel_pos, text = wm_star, fontsize = 18, align = (:left, :center))
## + ΔΘ
upper_ΔΘ = L"\left(S^{ΔΘ},\Theta^{*}+\Delta\Theta\right)"
text!(ax, S_ΔΘ[2]-0.005, ΔΘ_vec[2], text = upper_ΔΘ, fontsize = 18, align = (:right, :center))
lower_ΔΘ = L"\left(S^{-\Delta\Theta},\Theta^{*}-\Delta\Theta\right)"
text!(ax, S_ΔΘ[1]-0.005, ΔΘ_vec[1], text = lower_ΔΘ, fontsize = 18, align = (:right, :center))
## Stability
text!(ax, S[100], Θ[900], text = "Statically stable", color = density_grad[1])
text!(ax, S[900], Θ[100], text = "Statically unstable", color = density_grad[end],
      align = (:right, :center))
text!(ax, S[100], Θ[235], text = "Statically stable,\nunstable to cabbeling",
      color = density_grad[2])
arrows!(ax, [S[200]], [Θ[240]], [S[400]-S[200]], [Θ[100] - Θ[200]];
        lengthscale = 0.77, arrowcolor = density_grad[2], linecolor = density_grad[2])
arrows!(ax, [S[200]], [Θ[320]], [S[800]-S[200]], [Θ[700] - Θ[200]];
        lengthscale = 1.07, arrowcolor = density_grad[2], linecolor = density_grad[2])
## midpoint label for ΔΘ, if using rangebars
# text!(ax, water_parcel_pos[1], ΔΘ_vec[1] / 2, text = L"\Delta\Theta")
# text!(ax, Sₗ-0.005, ΔΘ_vec[2] / 2, text = L"\Delta\Theta", align = (:right, :center))
fig
