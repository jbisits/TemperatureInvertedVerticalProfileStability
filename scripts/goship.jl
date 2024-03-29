using .VerticalProfileStability
using MAT

## Density difference for all ocean profile data from GO-SHIP
oceans = ("atlantic", "indian", "pacific", "southern")
goship_analysis = joinpath(@__DIR__, "..", "data", "analysis", "goship.jld2")
ΔΘ_thres = (0.5, 1.0, 2.0, 3.0)

for ocean ∈ oceans

    files = glob("*.mat", joinpath(GOSHIP_DATADIR, ocean))
    @info "$(ocean) ocean"
    for ΔΘ ∈ ΔΘ_thres

        @info "- ΔΘ = $(ΔΘ)"
        res = goship_max_Δρ(files, ΔΘ)

        jldopen(goship_analysis, "a+") do file
            file["$(ΔΘ)/" * ocean] = res
        end

    end

end

## Create `all` entry that has all the `oceans` in one dictionary entry
goship_joined = joinpath(@__DIR__, "..", "data", "analysis", "goship_joined.jld2")
goship_output = jldopen(goship_analysis)

output_keys = ("Δρˢ", "Δρᶜ", "Θᵤ", "Θₗ", "pᵤ", "pₗ", "Sᵤ", "Sₗ", "lons", "lats")
for key ∈ keys(goship_output)

    temp_dict = Dict{String, Vector}(key => Vector{Union{Float64, Missing}}(undef, 0)
                                     for key ∈ output_keys)
    for ocean ∈ oceans

        for (i, data_key) ∈ enumerate(output_keys)

            append!(temp_dict[data_key], goship_output[key][ocean][data_key])

        end

    end

    jldopen(goship_joined, "a+") do file
        file[key] = temp_dict
    end
end

close(goship_output)
## Test function ouptut
goship_files = glob("*.mat", joinpath(GOSHIP_DATADIR, oceans[1]))
test_goship = goship_max_Δρ(goship_files[3:4], 2.0)

test_temp = collect(skipmissing(test_goship["Θₗ"]))
test_dd = collect(skipmissing(test_goship["Δρˢ"]))
scatter(test_temp, test_dd)
rm_extreme = findall(test_dd  .≤ 20)
fig, ax = scatter(test_temp[rm_extreme], test_dd[rm_extreme])
ylims!(ax, -0.1, 0.01)
fig

## Little data expoloration
goship_files = glob("*.mat", joinpath(GOSHIP_DATADIR, oceans[1]))
vars = matread(goship_files[4])["D_reported"]
keys(vars)
size(vars["lonlist"])
lons = vars["lonlist"]
lats = vars["latlist"]

scatter(lons, lats)

vec(vars["CTDprs"])[1]
## Measurements I want
obs = ("CTDSA", "CTDCT", "CTDprs", "lonlist", "latlist")

Sₐ, Θ, p = vars[obs[1]], vars[obs[2]], vars[obs[3]]
typeof(p)

p_vec_mat = (p)
p_vec_mat[1]
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
