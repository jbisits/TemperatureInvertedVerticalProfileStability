using .VerticalProfileStability
using MAT

ocean = ("atlantic", "indian", "pacific", "southern")
goship_files = glob("*.mat", joinpath(GOSHIP_DATADIR, ocean[1]))

ΔΘ_thres = (0.5, 1.0, 2.0, 3.0)

Main.VerticalProfileStability.MaximumDensityDifference.goship_max_Δρ(goship_files, ΔΘ_thres[3])

## Test function ouptut
test_goship = Main.VerticalProfileStability.MaximumDensityDifference.goship_max_Δρ(goship_files[1:2], 2.0)

test_temp = collect(skipmissing(test_goship["Θₗ"]))
test_dd = collect(skipmissing(test_goship["Δρˢ"]))
scatter(test_temp, test_dd)
rm_extreme = findall(test_dd  .≤ 20)
fig, ax = scatter(test_temp[rm_extreme], test_dd[rm_extreme])
ylims!(ax, -0.1, 0.01)
fig

## Little data expoloration
vars = matread(goship_files[1])["D_reported"]
keys(vars)
size(vars["lonlist"])
lons = vec(vars["lonlist"][2])
lats = vec(vars["latlist"])[1]

scatter(lons, lats)

vec(vars["CTDprs"])[1]
## Measurements I want
obs = ("CTDSA", "CTDCT", "CTDprs", "lonlist", "latlist")

Sₐ, Θ, p = vars[obs[1]][1][:, 1], vars[obs[2]][1][:, 1], vars[obs[3]][1][:, 1]

lines(Sₐ, Θ)
lines(gsw_rho.(Sₐ, Θ, p), p)
gsw_rho.(Sₐ, Θ, p)

##
vars[obs[1]][1]
findall(p .≤ 1000)
p[180]

## test `Dict` formation
test_keys = ("Δρˢ", "Δρᶜ", "Θᵤ", "Θₗ", "pᵤ", "pₗ", "Sᵤ", "Sₗ",
                                    "lons", "lats")
test_dict = Dict{String, Vector}(key => Float64[] for key ∈ test_keys)

test_keys[1]
