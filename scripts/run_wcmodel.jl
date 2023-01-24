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
Δt = 5 # 5 minute time-step, should remain stable

run_model(; Tᵤ, Tₗ, Sᵤ, Sᵣ, Sₗ, Sₘ, ref_density, savepath, Δt)
