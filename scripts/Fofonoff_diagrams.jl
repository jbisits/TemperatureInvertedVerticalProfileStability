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
# Save vals for other initial conditions
S_vals = Matrix{Float64}(undef, 11, length(T_increments))
T_vals = Matrix{Float64}(undef, 2, length(T_increments))
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
    Tₗ_2_find = findfirst(T .>= Tₗ_2)-1
    Sₗ_2_find = findfirst(ρ[:, Tₗ_2_find] .>= lower_isopycnal)-1
    Sₗ_2 = S[Sₗ_2_find]

    αₗ_2 = gsw_alpha(Sₗ_2, Tₗ_2, p_ref)
    βₗ_2 = gsw_beta(Sₗ_2, Tₗ_2, p_ref)
    m_2 = βₗ_2 / αₗ_2

    Tᵤ_2 = fill(Tᵤ_val + T_increment, 10)
    Tᵤ_2_idx = findfirst(T .> Tᵤ_2[1])
    Sᵤ_lower_isopycnal_idx = findfirst(ρ[:, Tᵤ_2_idx] .>= lower_isopycnal)
    Sᵤ_2_max = S[Sᵤ_lower_isopycnal_idx] - S_ic_dist_to_isopycnal
    s_range = range(Sᵤ_2_max - 0.07, Sᵤ_2_max, length = 10)

    T_vals[:, i] = [Tᵤ_2[1], Tₗ_2]
    S_vals[:, i] = vcat(s_range, Sₗ_2)

    inc = 0.01
    S_ = range(s_range[1] - inc, Sₗ_2 + inc, length = 1000)
    S_grid_ = S_ .* ones(length(S_))'
    T_ = range(Tᵤ_2[1] - inc, Tₗ_2 + inc, length = 1000)
    T_grid_ = T_' .* ones(length(T_))
    ρ_ = gsw_rho.(S_grid_, T_grid_, p_ref)

    tang_start = findfirst(S_ .>= s_range[4])
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
    scatter!(ax[i + 1], s_range, Tᵤ_2; color = ic_colour, markersize = 4)
end
fig
##
#save(joinpath(plotdir, "Fof_diagram_multiple_dwps.png"), fig)
##

## Threshold Δρ for given ΔΘ

## Isopycnal and linearised density about some deep water parcel
T = range(-2, 2, 200)
T_grid = T' .* ones(length(T))
S = range(34.4, 35, 200)
S_grid = S .* ones(length(S))'

p_ref = 0.0 #reference pressure
freezing_pt = gsw_ct_freezing.(S, p_ref, 1)

ρ = gsw_rho.(S_grid, T_grid, p_ref)

Sₗ, Tₗ = 34.7,  0.5
lower_isopycnal = gsw_rho(Sₗ, Tₗ, p_ref)

fig = Figure(resolution = (500, 500))
ax = Axis(fig[1, 1],
        xlabel = "Absolute salinity (g/kg)",
        ylabel = "Conservative temperature (∘C)")

contour!(ax, S, T, ρ,
        color = :black,
        linewidth = 2,
        levels = [lower_isopycnal],
        label = "Isopycnal (1027.71)")

αₗ = gsw_alpha(Sₗ, Tₗ, p_ref)
βₗ = gsw_beta(Sₗ, Tₗ, p_ref)
m = βₗ / αₗ

S_linear = range(34.51, Sₗ+0.05; length = 50)
Θ_linear = @. Tₗ + m * (S_linear - Sₗ)
scatter!(ax, [Sₗ], [Tₗ]; color = :red)
lines!(ax, S_linear, Θ_linear; color = :red)
fig

## Threshold Δρ for different values. Need to sort this out as its not right at the moment.

δΘ = [0.25, 0.5, 1, 2]
δS = @. Sₗ - (αₗ / βₗ) * (δΘ)
scatter!(ax, δS, Tₗ .- δΘ; color = :orange)
fig
save(joinpath(plotdir, "Δρ_thres_ex.png"), fig)

Tₗ_2 = Tₗ + 2
δΘ = 0:-0.01:-2
δS = @. Sₗ + (αₗ / βₗ) * (δΘ - Tₗ)
Δρ_thres = @. gsw_rho(δS, δΘ - Tₗ, p_ref) - gsw_rho(Sₗ, Tₗ, p_ref)

lines(δΘ, Δρ_thres)

## Looking at single δΘ for mulitple lower level temperatures

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

fig = Figure()
ax = Axis(fig[1, 1];
          xlabel = "Absolute salinity (g/kg)",
          ylabel = "Conservative temperature (ᵒC)")
contour!(ax, S, T, ρ;
         color = :black,
         linewidth = 2,
         levels = [lower_isopycnal],
         label = "Isopycnal (1027.71)")

#Θ = range(0, 20, length = 100)
Θ = LinRange(0.5, 20, 100)

finds_T = [findfirst(T .>= Θ_) for Θ_ ∈ Θ]
finds_S = [findfirst(ρ[:, finds_T_] .>= lower_isopycnal) for finds_T_ ∈ finds_T]
Sₐ = S[finds_S]
for i ∈ eachindex(Θ)
    scatter!(ax, [Sₐ[i]], [Θ[i]]; color = :red)
end
fig

## Δρ_thres
p_ref = 0.0
δΘ = [0.25, 0.5, 1.0, 2.0]
Δρ_thres_lower = Array{Float64}(undef, length(Θ), length(δΘ))
Δρ_thres_upper = similar(Δρ_thres_lower)

for (j, δΘ_) ∈ enumerate(δΘ)
    for i ∈ eachindex(Θ)

        αₗ = gsw_alpha(Sₐ[i], Θ[i], p_ref)
        βₗ = gsw_beta(Sₐ[i], Θ[i], p_ref)

        δS_lower = Sₐ[i] - (αₗ / βₗ) * (δΘ_)
        δS_upper = Sₐ[i] + (αₗ / βₗ) * (δΘ_)
        # scatter!(ax, [δS_lower], [Θ[i] - δΘ_]; color  = :orange)
        # scatter!(ax, [δS_upper], [Θ[i] + δΘ_]; color  = :green)
        Δρ_thres_lower[i, j] = gsw_rho(δS_lower, Θ[i] - δΘ_, p_ref) -
                               gsw_rho(Sₐ[i], Θ[i], p_ref)
        Δρ_thres_upper[i, j] = gsw_rho(δS_upper, Θ[i] + δΘ_, p_ref) -
                               gsw_rho(Sₐ[i], Θ[i], p_ref)

    end
end
fig

fig2 = Figure()
ax = Axis(fig2[1, 1];
          xlabel = "Θ (ᵒC) of lower level",
          xaxisposition = :top,
          ylabel = "Δρ_thres (kgm⁻³)",
          title = "Δρ threshold for ΔΘ threshold = $δΘ ᵒC against temperature of lower level"
          )
#lines!(ax, Θ, Δρ_thres_lower; label = "Upper water parcel colder")
#series!(ax, Θ, Δρ_thres_upper'; labels = string.(δΘ) .* "ᵒC")
series!(ax, Θ, Δρ_thres_lower'; labels = string.(δΘ) .* "ᵒC")
axislegend(ax; position = :rb)
fig2
#save(joinpath(plotdir, "Δρ_thres_multiple_deg.png"), fig2)
