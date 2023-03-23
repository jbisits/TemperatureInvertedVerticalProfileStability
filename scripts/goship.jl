using .VerticalProfileStability
using MAT

## Density difference for all ocean profile data
oceans = ("atlantic", "indian", "pacific", "southern")
goship_analysis = joinpath(@__DIR__, "..", "data", "analysis", "goship.jld2")
ΔΘ_thres = (0.5, 1.0, 2.0, 3.0)

for ocean ∈ oceans

    files = glob("*.mat", joinpath(GOSHIP_DATADIR, ocean))
    @info "$(ocean) ocean"
    for ΔΘ ∈ ΔΘ_thres

        res = goship_max_Δρ(files, ΔΘ)

        jldopen(goship_analysis, "a+") do file
            file["$(ΔΘ)/" * ocean] = res
        end

    end

end

## Test function ouptut
goship_files = glob("*.mat", joinpath(GOSHIP_DATADIR, oceans[1]))
test_goship = goship_max_Δρ(goship_files[1:2], 2.0)

test_temp = collect(skipmissing(test_goship["Θₗ"]))
test_dd = collect(skipmissing(test_goship["Δρˢ"]))
scatter(test_temp, test_dd)
rm_extreme = findall(test_dd  .≤ 20)
fig, ax = scatter(test_temp[rm_extreme], test_dd[rm_extreme])
ylims!(ax, -0.1, 0.01)
fig

## Little data expoloration
goship_files = glob("*.mat", joinpath(GOSHIP_DATADIR, oceans[1]))
vars = matread(goship_files[2])["D_reported"]
keys(vars)
size(vars["lonlist"])
lons = vec(vars["lonlist"])
lats = vec(vars["latlist"])

scatter(lons, lats)

vec(vars["CTDprs"])[1]
## Measurements I want
obs = ("CTDSA", "CTDCT", "CTDprs", "lonlist", "latlist")

Sₐ, Θ, p = vars[obs[1]], vars[obs[2]], vars[obs[3]]
typeof(p)
[p]
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
