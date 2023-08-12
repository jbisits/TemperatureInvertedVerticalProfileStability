using GibbsSeaWater, CairoMakie, SpecialFunctions

res = 1000
z = range(-50, 0; length = res)
offset = sum(extrema(z)) / 2
S₀ˡ, Θ₀ˡ = 34.7, 0.5
#S₀ᵘ, Θ₀ᵘ = 34.58565, -1.5
S₀ᵘ, Θ₀ᵘ = 34.564, -1.5
ΔS₀, ΔΘ₀ = S₀ᵘ -S₀ˡ, Θ₀ᵘ - Θ₀ˡ

"""
    initial_tracer_heaviside(z::AbstractVector, C::Number, ΔC::Number)
Modified Heaviside function for initial condition of a tracer with depth. The offset of the
Heaviside function is calculated from the `extrema` of the depth array `z`.

Function arguments:

- `z` depth of the domain;
- `C` tracer value in deeper part of the step;
- `ΔC` difference in tracer between the steps.
"""
function initial_tracer_heaviside(z::AbstractVector, C::Number, ΔC::Number)

    initial_values = similar(z)
    offset = sum(extrema(z)) / 2
    for i ∈ eachindex(z)
        initial_values[i] = z[i] - offset < 0 ? C : C + ΔC
    end

    return initial_values

end

"""
    function tracer_solution(z::AbstractVector, C::Number, ΔC::Number,
                             t::Union{AbstractVector, Number})
Solution to the heat equation for a tracer concentration field `C` subject to initial
conditions that are a Heaviside (or modified Heaviside) step function. Returned is a matrix
with the solutions at each depth level down the column and the each column being a time step.

Function arguments:

- `z` depth of the domain;
- `C` tracer value in deeper part of the step;
- `ΔC` difference in tracer between the steps;
- `κ` the diffusivity of the tracer;
- `t` time at which to evaluate the solution, can either be a `Vector` of time steps or just
a number.
"""
function tracer_solution(z::AbstractVector, C::Number, ΔC::Number, κ::Number,
                         t::Union{AbstractVector, Number})

    offset = sum(extrema(z)) / 2
    res = typeof(t) <: Number ? Array{Float64}(undef, length(z)) :
                                Array{Float64}(undef, length(z), length(t))

    if typeof(t) <: Number
        for i ∈ eachindex(z)
            res[i] = C + 0.5 * ΔC * (1 + erf((z[i] - offset) / sqrt(4 * κ * t)))
        end
    else
        for i ∈ eachindex(z), j ∈ eachindex(t)
            res[i, j] = C + 0.5 * ΔC * (1 + erf((z[i] - offset) / sqrt(4 * κ * t[j])))
        end
    end

    return res

end

##
t = 1:3:10
κₜ, κₛ = 1e-2, 1e-2
initial_Θ = initial_tracer_heaviside(z, Θ₀ˡ, ΔΘ₀)
Θ_sol = tracer_solution(z, Θ₀ˡ, ΔΘ₀, κₜ, t)
initial_S = initial_tracer_heaviside(z, S₀ˡ, ΔS₀)
S_sol = tracer_solution(z, S₀ˡ, ΔS₀, κₛ, t)
initial_σ₀ = gsw_sigma0.(initial_S, initial_Θ)
σ₀_sol = gsw_sigma0.(S_sol, Θ_sol)
lines(Θ_sol[:, 4], z)

##
t = 20
κₜ, κₛ = 1e-2, 1e-2
initial_Θ = initial_tracer_heaviside(z, Θ₀ˡ, ΔΘ₀)
Θ_sol = tracer_solution(z, Θ₀ˡ, ΔΘ₀, κₜ, t)
initial_S = initial_tracer_heaviside(z, S₀ˡ, ΔS₀)
S_sol = tracer_solution(z, S₀ˡ, ΔS₀, κₛ, t)
initial_σ₀ = gsw_sigma0.(initial_S, initial_Θ)
σ₀_sol = gsw_sigma0.(S_sol, Θ_sol)

fontsize = 22
labelsize = 16
fig = Figure(size = (900, 400); fontsize)
ax = [Axis(fig[1, i]) for i ∈ 1:2]
lines!(ax[1], initial_S, z; color = (:blue, 0.5), label = "Initial salinity")
lines!(ax[1], S_sol, z; color = :blue, linestyle = :dash, label = "Salinity after mixing")
xlims!(ax[1], 34.55, 34.72)
ax[1].xlabel = "Salinity (gkg⁻¹)"
ax[1].xlabelcolor = :blue
ax[1].xticklabelcolor = :blue
ax[1].ylabel = "z (m)"
ax2 = Axis(fig[1, 1];
           xaxisposition = :top,
           xticklabelcolor = :red,
           xlabel = "Θ (°C)",
           xlabelcolor = :red,
           title = "Temperature and salinity profile")
lines!(ax2, initial_Θ, z; color = (:red, 0.5), label = "Initial temeperature")
lines!(ax2, Θ_sol, z; color = :red, linestyle = :dash, label = "Temperature after mixing")
lines!(ax[2], initial_σ₀, z; color = (:black, 0.5), label = "Initial density")
lines!(ax[2], σ₀_sol, z; color = :black, linestyle = :dash, label = "Density after mixing")
ax[2].title = "Density profile"
ax[2].xlabel = "σ₀ (kgm⁻³)"
ax[2].xticklabelrotation = π/4
xlims!(ax[2], minimum(initial_σ₀)-0.001, maximum(σ₀_sol)+0.001)
axislegend(ax2; labelsize)
axislegend(ax[1], position = :lb; labelsize)
axislegend(ax[2]; labelsize)
hideydecorations!(ax[2], grid = false)

linkyaxes!(ax[1], ax[2])
colsize!(fig.layout, 1, Relative(3/5))
fig

##
save("plots/profile2.png", fig)
