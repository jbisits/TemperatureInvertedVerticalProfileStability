using CairoMakie, ColorSchemes, GibbsSeaWater

## Setting up initial conditions for the model using a Fofonoff diagram.
T = range(-2, 22, 5000)
T_grid = T' .* ones(length(T))
S = range(34.4, 40, 5000)
S_grid = S .* ones(length(S))'

p_ref = 0 #reference pressure
freezing_pt = gsw_ct_freezing.(S, p_ref, 1)

ρ = gsw_rho.(S_grid, T_grid, p_ref)

Sₗ, Tₗ = 34.7,  0.5
Tᵤ_val = -1.85
lower_isopycnal = gsw_rho(Sₗ, Tₗ, p_ref)
ic_colour = reverse(get(ColorSchemes.viridis, range(0, 1, length = 10)))

T_increments = 5 .* [0, 1, 2, 3, 4]
##
fig = Figure(resolution = (1400, 1400))
titles = reshape(["Θ-S diagram with multiple deep water masses"
                 ["Θ-S diagram\ndeep water conservative temperature = $(Tₗ + T_increment)ᵒC"
                 for T_increment ∈ T_increments]], (2, 3))
ax = [Axis(fig[i, j],
        xlabel = "Absolute salinity (g/kg)",
        ylabel = "Conservative temperature (∘C)",
        title = titles[i, j]) for i ∈ 1:2, j ∈ 1:3]

contour!(ax[1], S, T, ρ,
        color = :black,
        linewidth = 2,
        levels = [lower_isopycnal],
        label = "Isopycnal (1027.71)")
lines!(ax[1], S, freezing_pt, linestyle = :dash, color = :blue, label = "Freezing point")
axislegend(ax[1], position = :lt)
fig
# Distance that the initial condition with highest salinity is from the `lower_isopycnal`.
upper_level_T_idx = findfirst(T .> Tᵤ_val)
S_max_sal_idx = findfirst(ρ[:, upper_level_T_idx] .>= lower_isopycnal)
S_ic_dist_to_isopycnal = S[S_max_sal_idx] - 34.57

# Increment and plot on separate panels
deep_water_colour = [:red, :orange, :blue, :green, :magenta]
for (i, T_increment) ∈ enumerate(T_increments)

    Tₗ_2 = Tₗ + T_increment
    Tₗ_2_find = findfirst(T .>= Tₗ_2)
    Sₗ_2_find = findfirst(ρ[:, Tₗ_2_find] .>= lower_isopycnal)
    Sₗ_2 = S[Sₗ_2_find]

    αₗ_2 = gsw_alpha(Sₗ_2, Tₗ_2, p_ref)
    βₗ_2 = gsw_beta(Sₗ_2, Tₗ_2, p_ref)
    m_2 = βₗ_2 / αₗ_2

    Tᵤ_2 = fill(Tᵤ_val + T_increment, 10)
    Tᵤ_2_idx = findfirst(T .> Tᵤ_2[1])
    Sᵤ_lower_isopycnal_idx = findfirst(ρ[:, Tᵤ_2_idx] .>= lower_isopycnal)
    Sᵤ_2_max = S[Sᵤ_lower_isopycnal_idx] - S_ic_dist_to_isopycnal
    s_range = range(Sᵤ_2_max - 0.1, Sᵤ_2_max, length = 10)

    inc = 0.01
    S_ = range(s_range[1] - inc, Sₗ_2 + inc, length = 1000)
    S_grid_ = S_ .* ones(length(S_))'
    T_ = range(Tᵤ_2[1] - inc, Tₗ_2 + inc, length = 1000)
    T_grid_ = T_' .* ones(length(T_))
    ρ_ = gsw_rho.(S_grid_, T_grid_, p_ref)

    tang_start = findfirst(S_ .>= s_range[5])
    tang_length = tang_start:findfirst(S_ .>= Sₗ_2)
    S_tangent = S_[tang_length]
    tangent = @. Tₗ_2 + m_2 * (S_tangent - Sₗ_2)

    scatter!(ax[1], [Sₗ_2], [Tₗ_2]; color = deep_water_colour[i],
            label = "Deep water mass")

    contour!(ax[i + 1], S_, T_, ρ_,
            color = :black,
            linewidth = 2,
            levels = [lower_isopycnal],
            label = "Isopycnal (1027.71)")
    scatter!(ax[i + 1], [Sₗ_2], [Tₗ_2]; color = deep_water_colour[i],
            label = "Deep water mass")
    lines!(ax[i + 1], S_tangent, tangent; color = deep_water_colour[i],
            label = "Linearised density about\ndeep water mass")
    scatter!(ax[i + 1], s_range, Tᵤ_2; color = ic_colour, markersize = 6)
end
fig
##
save(joinpath(plotdir, "Fof_diagram_multiple_dwps.png"), fig)
