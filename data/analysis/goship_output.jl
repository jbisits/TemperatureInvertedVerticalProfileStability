using .VerticalProfileStability

const GOSHIP_DATA = joinpath(@__DIR__, "goship.jld2")
const GOSHIP_JOINED = joinpath(@__DIR__, "goship_joined.jld2")

gd = jldopen(GOSHIP_DATA)

ΔΘ_vals = (0.5, 1.0, 2.0, 3.0)
oceans = ("atlantic", "indian", "pacific", "southern")

fig = Figure()
ax = Axis(fig[1, 1];
          xlabel = "Longitude",
          ylabel = "Latitude")

for key ∈ keys(gd["1.0"])

    Θᵤ = collect(skipmissing(gd["1.0"][key]["Θᵤ"]))
    Θₗ = collect(skipmissing(gd["1.0"][key]["Θₗ"]))
    find_inv = findall(Θᵤ .< Θₗ)
    lons = gd["1.0"][key]["lons"][find_inv]
    lats = gd["1.0"][key]["lats"][find_inv]
    scatter!(ax, lons, lats; label = key)

end
Legend(fig[2, 1], ax, orientation = :horizontal)
fig

fig2 = Figure()
ax2 = Axis(fig2[1, 1];
          xlabel = "Θₗ",
          ylabel = "Δρ_static")
for ΔΘ_key ∈ keys(gd)
    for key ∈ keys(gd["1.0"])
        Δρˢ = collect(skipmissing(gd[ΔΘ_key][key]["Δρˢ"]))
        Θᵤ = collect(skipmissing(gd[ΔΘ_key][key]["Θᵤ"]))
        Θₗ = collect(skipmissing(gd[ΔΘ_key][key]["Θₗ"]))
        find_inv = Θᵤ .< Θₗ
        scatter!(ax2, Θₗ[find_inv], Δρˢ[find_inv]; label = key)
    end
end
Legend(fig2[2, 1], ax, orientation = :horizontal)
fig2
ylims!(ax2, -0.1, 0.01)
xlims!(ax2, -1.8, 15)
fig2
##
close(gd)

## Joined data
gdj = jldopen(GOSHIP_JOINED)
gdj["0.5"]["all"]


##
close(gdj)
