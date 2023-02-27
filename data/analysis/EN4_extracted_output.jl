using .VerticalProfileStability

const EN4_OUTPUT = joinpath(@__DIR__, "EN4_extracted.jld2")

en4_data = jldopen(EN4_OUTPUT)
data_keys = keys(en4_data)
dict_keys = keys(en4_data["EN4_2007"])
merged_dict = Dict{String, Any}()
for key ∈ dict_keys

    dummy = en4_data["EN4_2007"][key]
    for year ∈ data_keys
        dummy = vcat(dummy, en4_data[year][key])
    end
    push!(merged_dict, key => dummy)
end
merged_dict

Θᵤ = merged_dict["Θᵤ"]
Sᵤ = merged_dict["Sᵤ"]
pᵤ = merged_dict["pᵤ"]
Θₗ = merged_dict["Θₗ"]
Sₗ = merged_dict["Sₗ"]
pₗ = merged_dict["pₗ"]
Δρˢ = merged_dict["Δρˢ"]
find_inversion = findall(Θᵤ .< Θₗ) #∩ findall(Θₗ .≤ 10)

fig = Figure(size = (500, 500))
ax = Axis(fig[1, 1];
          xlabel = "Θₗ (°C)",
          ylabel = "Δρ (kgm⁻³)")
scatter!(ax, Θₗ[find_inversion], Δρˢ[find_inversion])
Θ_range = range(-2, 10; length = 100)
# Sₗ_mean = mean(Sₗ[find_inversion])
# p̄_mean = mean(0.5 .* (pₗ[find_inversion] .+ pᵤ[find_inversion]))
Sₗ_mean = 35
p̄_mean = 300
αₗ = gsw_alpha.(Sₗ_mean, Θ_range, p̄_mean)
βₗ = gsw_beta.(Sₗ_mean, Θ_range, p̄_mean)

Δρ_thres = @. gsw_rho(Sₗ_mean - (αₗ / βₗ) * 1, Θ_range - 1, p̄_mean) - gsw_rho(Sₗ_mean, Θ_range, p̄_mean)
lines!(ax, Θ_range, Δρ_thres)
fig
