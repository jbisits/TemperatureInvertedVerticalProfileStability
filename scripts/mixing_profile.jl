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
    lines!(a, plot_vars[i], z)
    if isodd(i)
        a.xticklabelcolor = "blue"
        lines!(ax2[counter], plot_vars2[counter], z; color = (:red, 0.5))
        ax2[counter].xticklabelcolor = "red"
        ax2[counter].xaxisposition = :top
        counter += 1
    end
    #a.title = titles[i]
end
linkxaxes!()
linkxaxes!(ax[2], ax[4])
fig
