using GibbsSeaWater, CairoMakie

fig = Figure(size = (1200, 1200))
ax = [Axis(fig[i, j]) for j ∈ 1:2, i ∈ 1:2]
ax2 = [Axis(fig[i, 1]) for i ∈ 1:2]
titles = ("Initial salt and temperature profiles", "Initial density profile",
          "Post mixing salt and temperature profiles", "Post mixing density profile")
res = 1000
z = range(-50, 0; length = res)
z_mid = sum(extrema(z)) / 2
Sᵤ, Sₗ = 34.58565, 34.7
ΔS = (Sᵤ - Sₗ) / 2
S = vcat(fill(Sₗ, round(Int, res / 2)), fill(Sᵤ, round(Int, res / 2)))
S_mixed = @. ΔS * tanh(z - z_mid) + Sₗ + ΔS
Θᵤ, Θₗ = -1.5, 0.5
ΔΘ = (Θᵤ - Θₗ) / 2
Θ = vcat(fill(Θₗ, round(Int, res / 2)), fill(Θᵤ, round(Int, res / 2)))
Θ_mixed = @. ΔΘ * tanh(z - z_mid) + Θₗ + ΔΘ
σ₀ = round.(gsw_sigma0.(S, Θ); digits = 4)
lines(σ₀, z)
lines(S, z)
lines(S_mixed, z)
lines(Θ_mixed, z)

σ₀_mixed = round.(gsw_sigma0.(S_mixed, Θ_mixed); digits = 4)
lines(σ₀_mixed, z)

plot_vars = (S, σ₀, S_mixed, σ₀_mixed)
plot_vars2 = (Θ, Θ_mixed)
counter = 1
for (i, a) ∈ enumerate(ax)
    lines!(a, plot_vars[i], z;
          color = isodd(i) ? :blue : :orange,
          linestyle = isodd(i) ? :dot : :solid)
    a.xlabel = "Salinity (gkg⁻¹)"
    a.ylabel = "z (m)"
    if isodd(i)
        a.xticklabelcolor = "blue"
        lines!(ax2[counter], plot_vars2[counter], z; color = (:red, 0.5), linestyle = :dash)
        ax2[counter].xticklabelcolor = "red"
        ax2[counter].xaxisposition = :top
        ax2[counter].xlabel = "Θ (°C)"
        ax2[counter].xlabelcolor = :red
        ax2[counter].title = titles[i]
        xlims!(a, 34.55, 34.72)
        a.xlabelcolor = :blue
        counter += 1
    end
    if iseven(i)
        a.title = titles[i]
        a.xlabel = "σ₀ (kgm⁻³)"
        xlims!(a, 27.705, 27.7135)
        a.xticklabelrotation = π/4
        hideydecorations!(a, grid = false)
    end
end
linkxaxes!(ax[1], ax[3])
linkxaxes!(ax[2], ax[4])
colsize!(fig.layout, 1, Relative(3/5))
fig

##
using SpecialFunctions
heaviside(z) = z ≤ 0 ? 0 : 1
heaviside(z, upper, lower) = z ≤ 0 ? lower : upper
heaviside(z, z_mid, upper, lower) = z ≤ z_mid ? lower : upper
Θ_evolution(z, t, ΔΘ, Θ_offset; κ = 1) = 0.5 * ΔΘ * (1 + erf(z / sqrt(4 * κ * t))) + Θ_offset
S_evolution(z, t, ΔS, S_offset; κ = 1) = 0.5 * ΔS * (1 + erf(z / sqrt(4 * κ * t))) + S_offset
z = range(-1, 1; length = 100)
z2 = range(-2, 0; length = 100)
t = range(0.1, 1; length = 5)
S_res = Array{Float64}(undef, length(z), length(t))
Θ_res = Array{Float64}(undef, length(z), length(t))
for (i, z_) ∈ enumerate(z2), (j, t_) ∈ enumerate(t)
    S_res[i, j] = S_evolution(z_ + 1, t_, 0.11435, 34.7, κ = 1e-2)
    Θ_res[i, j] = Θ_evolution(z_ + 1, t_, -2, 0.5, κ = 1e-2)
end
S_res = hcat(heaviside.(z2, -1, 34.58565, 34.7))
Θ_res = hcat(heaviside.(z2, -1, -1.5, 0.5), Θ_res)
fig, ax = lines(Θ_res[:, 1], z2)
for sol ∈ eachcol(Θ_res)
    lines!(sol, z2)
end
fig
