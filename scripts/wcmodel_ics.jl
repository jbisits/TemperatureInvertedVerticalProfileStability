## Generate more water column model initial conditions for higher temperatures. To do this
#  I increase the overall temperature in each layer by the same amount (e.g. 1ᵒC) then find
#  salinity initial conditions that are the same distance from the isopycnal as what the
#  model has been run with.

using CairoMakie, ColorSchemes, GibbsSeaWater

## Setting up initial conditions for the model using a Fofonoff diagram.
T = range(-2, 9, 5000)
T_grid = T' .* ones(length(T))
S = range(34.4, 36, 5000)
S_grid = S .* ones(length(S))'

p_ref = 0 #reference pressure
freezing_pt = gsw_ct_freezing.(S, p_ref, 1)

ρ = gsw_rho.(S_grid, T_grid, p_ref)

Sₗ, Tₗ = 34.7,  0.5
lower_isopycnal = gsw_rho(Sₗ, Tₗ, p_ref)

αₗ = gsw_alpha(Sₗ, Tₗ, p_ref)
βₗ = gsw_beta(Sₗ, Tₗ, p_ref)
m = βₗ / αₗ

tang_start = 300
tang_length = tang_start:findfirst(S .> Sₗ)
S_tangent = S[tang_length]
tangent = @. Tₗ + m * (S_tangent - Sₗ)

Tᵤ_val = -1.85
Tᵤ = fill(Tᵤ_val, 10)
s_range = range(34.47, 34.57, length = 10)
ic_colour = reverse(get(ColorSchemes.viridis, range(0, 1, length = 10)))

ic_plot = Figure(resolution = (1000, 1000))

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
#lines!(ax, Sₗ:0.05:S[end], 0.5 * ones(length(Sₗ:0.05:S[end]));
#         color = :red, linestyle = :dot, label = "Salinity in lower layer")
scatter!(ax, s_range, Tᵤ, color = ic_colour)
#axislegend(ax, position = :lt)
Legend(ic_plot[1, 2], ax)
ic_plot

# Distance that the initial condition with highest salinity is from the `lower_isopycnal`.
upper_level_T_idx = findfirst(T .> Tᵤ_val)
S_max_sal_idx = findfirst(ρ[:, upper_level_T_idx] .>= lower_isopycnal)
S_ic_dist_to_isopycnal = S[S_max_sal_idx] - 34.57

# Create another 4 sets of initial conditions each increasing by 2ᵒC.
T_increments = [2, 4, 6, 8]
for T_increment ∈ T_increments
    Tₗ_2 = Tₗ + T_increment
    Tₗ_2_find = findfirst(T .>= Tₗ_2)
    Sₗ_2_find = findfirst(ρ[:, Tₗ_2_find] .>= lower_isopycnal)
    Sₗ_2 = S[Sₗ_2_find]

    αₗ_2 = gsw_alpha(Sₗ_2, Tₗ_2, p_ref)
    βₗ_2 = gsw_beta(Sₗ_2, Tₗ_2, p_ref)
    m_2 = βₗ_2 / αₗ_2

    tang_start += 300
    tang_length = tang_start:findfirst(S .>= Sₗ_2)
    S_tangent = S[tang_length]
    tangent = @. Tₗ_2 + m_2 * (S_tangent - Sₗ_2)

    Tᵤ_2 = fill(Tᵤ_val + T_increment, 10)
    Tᵤ_2_idx = findfirst(T .> Tᵤ_2[1])
    Sᵤ_lower_isopycnal_idx = findfirst(ρ[:, Tᵤ_2_idx] .>= lower_isopycnal)
    Sᵤ_2_max = S[Sᵤ_lower_isopycnal_idx] - S_ic_dist_to_isopycnal
    s_range_2 = range(Sᵤ_2_max - 0.1, Sᵤ_2_max, length = 10)

    scatter!(ax, [Sₗ_2], [Tₗ_2], color = :red, label = "Deep water mass")
    lines!(ax, S_tangent, tangent, color = :red,
            label = "Linearised density about\ndeep water mass")
    scatter!(ax, s_range_2, Tᵤ_2, color = ic_colour)
end
ic_plot

# Now can just write a function to find the values for `Tᵤ` and `Tₗ`,
# `Sᵤ`, `Sₗ`, `Sₘ` and `Sᵣ`. I think that should do everthing correctly. Although the model
# does not use density referenced to so might need adjustmetns when using the isopycnal
# for the initial conditions salinity distance?

## Zoom in on each part where an IC is

## Tₗ = 0.5
xlims!(ax, 34.4, 34.8)
ylims!(ax, -2, 1)
ic_plot

## Tₗ = 2.5
xlims!(ax, 34.51, 34.9)
ylims!(ax, 0, 3)
ic_plot

## Tₗ = 4.5
xlims!(ax, 34.7, 35.2)
ylims!(ax, 2, 5)
ic_plot

## Tₗ = 6.5
xlims!(ax, 34.9, 35.5)
ylims!(ax, 4, 7)
ic_plot

## Tₗ = 8.5
xlims!(ax, 35.3, 36)
ylims!(ax, 6, 9)
ic_plot
