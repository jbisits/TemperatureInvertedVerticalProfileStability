using CairoMakie, ColorSchemes, GibbsSeaWater

## Setting up initial conditions for the model using a Fofonoff diagram.
T = range(-2, 1, 200)
T_grid = T' .* ones(length(T))
S = range(34.4, 34.75, 200)
S_grid = S .* ones(length(S))'

p_ref = 0 #reference pressure
freezing_pt = gsw_ct_freezing.(S, p_ref, 1)

ρ = gsw_rho.(S_grid, T_grid, p_ref)

Sₗ, Tₗ = 34.7,  0.5
lower_isopycnal = gsw_rho(Sₗ, Tₗ, p_ref)

αₗ = gsw_alpha(Sₗ, Tₗ, p_ref)
βₗ = gsw_beta(Sₗ, Tₗ, p_ref)
m = βₗ / αₗ

tang_start = 56
tang_length = tang_start:findfirst(S .> Sₗ)
S_tangent = S[tang_length]
tangent = @. Tₗ + m * (S_tangent - Sₗ)

Tᵤ_val = -1.85
Tᵤ = fill(Tᵤ_val, 10)
s_range = range(34.47, 34.57, length = 10)
ic_colour = reverse(get(ColorSchemes.viridis, range(0, 1, length = 10)))

ic_plot = Figure(resolution = (500, 500))

ax = Axis(ic_plot[1, 1],
        xlabel = "Absolute salinity (g/kg)",
        ylabel = "Conservative temperature (∘C)",
        title = "Fofonoff diagram - used to find the\nsalinity initial conditions in mixed layer")

contour!(ax, S, T, ρ,
        color = :black,
        linewidth = 2,
        levels = [lower_isopycnal],
        label = "Isopycnal (1027.71)")
lines!(ax, S, freezing_pt, linestyle = :dash, color = :blue, label = "Freezing point")
scatter!(ax, [Sₗ], [Tₗ], color = :red, label = "Deep water mass")
lines!(ax, S_tangent, tangent, color = :red,
        label = "Linearised density about\ndeep water mass")
scatter!(ax, s_range, Tᵤ, color = ic_colour)
axislegend(ax, position = :lt)
ic_plot

## Setting up initial conditions for the model using a Fofonoff diagram.
T = range(19, 22, 200)
T_grid = T' .* ones(length(T))
S = range(36, 37, 200)
S_grid = S .* ones(length(S))'

p_ref = 0 #reference pressure

ρ = gsw_rho.(S_grid, T_grid, p_ref)

Sₗ, Tₗ = 36.8, 21
lower_isopycnal = gsw_rho(Sₗ, Tₗ, p_ref)

αₗ = gsw_alpha(Sₗ, Tₗ, p_ref)
βₗ = gsw_beta(Sₗ, Tₗ, p_ref)
m = βₗ / αₗ

tang_start = 15
tang_length = tang_start:findfirst(S .> Sₗ)
S_tangent = S[tang_length]
tangent = @. Tₗ + m * (S_tangent - Sₗ)

Tᵤ_val = 19.15
Tᵤ = fill(Tᵤ_val, 10)
s_range = range(36.03, 36.13, length = 10)
ic_colour = reverse(get(ColorSchemes.viridis, range(0, 1, length = 10)))

ic_plot = Figure(resolution = (900, 900))

ax = Axis(ic_plot[1, 1],
        xlabel = "Absolute salinity (g/kg)",
        ylabel = "Conservative temperature (∘C)",
        title = "Fofonoff diagram - used to find the\nsalinity initial conditions in mixed layer")

contour!(ax, S, T, ρ,
        color = :black,
        linewidth = 2,
        levels = [lower_isopycnal],
        label = "Isopycnal (1027.71)")
scatter!(ax, [Sₗ], [Tₗ], color = :red, label = "Deep water mass")
lines!(ax, S_tangent, tangent, color = :red,
        label = "Linearised density about\ndeep water mass")
scatter!(ax, s_range, Tᵤ, color = ic_colour)
axislegend(ax, position = :lt)
ic_plot
