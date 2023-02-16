using .VerticalProfileStability
using MAT, JLD2

const ARGO_DATA = joinpath(@__DIR__, "..", "data", "observations", "ARGO_JJASO",
                           "Argo_JJASO.mat")

## Quick look
file = matopen(ARGO_DATA)
varnames = keys(file)
vars = matread(ARGO_DATA)
close(file)

## I have not looked at too many but there does seem to be a few that will suffice. Not
#  sure what the variables are but likely potential temperature and practical salinity.

lats = round.(vec(vars["lat"]); digits = 2)
lon = vars["lon"]
S = vars["sal"]
T = vars["temp"]
p = vars["press"]

profile = 4000
lines(reshape(S[profile], :), reshape(T[profile], :))

## Extract maximum density differences
ΔΘ_thres = [[0.5, 1.0], [1.0, 2.0], [2.0, 3.0], 3.0]
Argo_ouput = joinpath(@__DIR__, "..", "data", "analysis", "ARGO_extracted.jld2")
for ΔΘ ∈ ΔΘ_thres

    @info "$(ΔΘ) threshold"
    argo = argo_max_Δρ(ARGO_DATA, ΔΘ)

    @info "Saving"
    jldopen(Argo_ouput, "a+") do file

        file["ΔΘ_thres_$(ΔΘ)"] = argo

    end

end
