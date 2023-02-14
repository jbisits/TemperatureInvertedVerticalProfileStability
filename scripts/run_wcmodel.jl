using .VerticalProfileStability

## Default parameters
run_model()

## Warmer temperature run. S_vals and T_vals come from `Fofonoff_diagrams.jl`

Tᵤ = T_vals[1, end]
Tₗ = T_vals[2, end]
Sᵤ = S_vals[1, end]
Sᵣ = S_vals[end-1, end]
Sₗ = S_vals[end, end]
Sₘ = S_vals[end, end] + 0.1
ref_density = gsw_rho(Sₗ, Tₗ, 100)
savepath = joinpath(sim_datadir, "initial_Θ_18_15_20_5")
Δt = 5

run_model(; Tᵤ, Tₗ, Sᵤ, Sᵣ, Sₗ, Sₘ, ref_density, savepath, Δt)

## Initial conditions based on Δρ_thres, I want three initial conditions which are either
# side of and on the critical salinity value found by S = Sₗ ± (αₗ / βₗ)ΔΘ.
p_ref = 0
S_increment = 0.002
num_ics = 3

## ΔΘ = 0.5
ΔΘ_vals = [0.5, 1.0, 2.0]
Sₗ, Tₗ = 34.7, 0.5
αₗ, βₗ = gsw_alpha(Sₗ, Tₗ, p_ref), gsw_beta(Sₗ, Tₗ, p_ref)

for ΔΘ ∈ ΔΘ_vals[1:2]

    Tᵤ = Tₗ - ΔΘ
    S_critical = Sₗ - (αₗ / βₗ) * ΔΘ
    Sᵤ = S_critical - S_increment
    Sᵣ = S_critical + S_increment
    savepath = joinpath(sim_datadir, "initial_ΔΘ_$(ΔΘ)")
    # mkdir(savepath)
    # Run Isopycnal and linearised density about some deep water parcel in
    # `Fofonoff_diagrams` to set up the axis and view the initial conditions
    #scatter!(ax, [Sᵤ, S_critical, Sᵣ], Tᵤ .* ones(3))
    run_model(; Sᵤ, Sᵣ, Sₗ, Tₗ, Tᵤ, savepath, num_ics)

end
