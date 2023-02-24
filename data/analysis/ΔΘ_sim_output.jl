using JLD2, ColorSchemes
using .VerticalProfileStability

############################################################################################
## Choose which simulation
############################################################################################
simulations = ("initial_ΔΘ_0.5", "initial_ΔΘ_1.0", "initial_ΔΘ_2.0")
sim_num = 3 # Change this to look at other simulations
sim_output = joinpath(SIM_DATADIR, simulations[sim_num])

############################################################################################
## Open and load data from chosen simulation, then extract initial conditions
############################################################################################
saved_ts = jldopen(joinpath(sim_output, "output_timeseries.jld2"))
t, S_ts, T_ts, κ_ts, Δρ_vals = saved_ts["t"], saved_ts["S_ts"], saved_ts["T_ts"],
                               saved_ts["κ_ts"], saved_ts["Δρ_ts"]
close(saved_ts)

z = -500:5:-5 # This can be accessed from the saved data but much faster this way

saved_params = jldopen(joinpath(sim_output, "model_params.jld2"))
Sᵤ, Sₗ, Sᵣ, Sₘ, Θᵤ, Θₗ, num_ics = saved_params["Sᵤ"], saved_params["Sₗ"], saved_params["Sᵣ"],
                                  saved_params["Sₘ"], saved_params["Tᵤ"], saved_params["Tₗ"],
                                  saved_params["num_ics"]
close(saved_params)
p_ref = 0

S₀, T₀ = [S_ts[:, 1, i] for i ∈ 1:length(S_ts[1, 1, :])],
         [T_ts[:, 1, i] for i ∈ 1:length(T_ts[1, 1, :])]
sal_ic_vals_string = string.(round.([S₀[i][end] for i ∈ 1:num_ics]; digits = 3))

############################################################################################
## Initial conditions for the simulation
############################################################################################

# Axis and constant setups
Θ_range = (-2, 1)
T = range(Θᵤ - 0.25, Θₗ + 0.25, 200)
T_grid = T' .* ones(length(T))
S = range(S₀[1][end] - 0.01, 34.71, 200)
S_grid = S .* ones(length(S))'

ρ = gsw_rho.(S_grid, T_grid, p_ref)

lower_isopycnal = gsw_rho(Sₗ, Θₗ, p_ref)

αₗ = gsw_alpha(Sₗ, Θₗ, p_ref)
βₗ = gsw_beta(Sₗ, Θₗ, p_ref)
m = βₗ / αₗ

tang_length = 7:findfirst(S .> Sₗ)
S_tangent = S[tang_length]
tangent = @. Θₗ + m * (S_tangent - Sₗ)

ic_colour = reverse(get(ColorSchemes.viridis, range(0, 1, length = num_ics)))

ic_colour = reverse(get(ColorSchemes.viridis, range(0, 1, length = num_ics * 3)))
ic_colour = reshape(ic_colour, 3, 3)
ΔΘ_colour = [:blue, :orange, :green]
ic_plot = Figure(resolution = (1000, 1000))
xlabs = ["Temperature (∘C)" "Salinity (g/kg)"; "Inital ΔΘ (°C) at interface " "Salinity (g/kg)"]
ylabs = ["Depth (m)" "Depth (m)"; "Δρ (kgm⁻³)" "Temperature (°C)"]
titles = ["(a) Temperature initial conditions" "(b) Salinity initial conditions";
          "(c) Density difference threshold " "(d) Fofonoff diagram of initial conditions"]
ax = [Axis(ic_plot[j, i],
        xlabel = xlabs[j, i],
        xaxisposition = (j == 1 ? :top : :bottom),
        ylabel = ylabs[j, i],
        title = titles[j, i]) for i ∈ 1:2, j ∈ 1:2]
linkyaxes!(ax[1, 1], ax[2, 1])
linkxaxes!(ax[2, 1], ax[2, 2])
# Fofonoff diagram
contour!(ax[2, 2], S, T, ρ,
        color = :black,
        linewidth = 2,
        levels = [lower_isopycnal],
        label = "Isopycnal (1027.71)")
lines!(ax[2, 2], S_tangent, tangent; color = :red,
       label = "Linearised density about\ndeep water mass")
scatter!(ax[2, 2], [Sₗ], [Θₗ], color = :red, label = "Deep water mass")

ic_plot
## Plotting
for (j, sim) ∈ enumerate(simulations)
    sim_output_ = joinpath(SIM_DATADIR, sim)
    saved_ts_ = jldopen(joinpath(sim_output_, "output_timeseries.jld2"))
    S_ts_, T_ts_ = saved_ts_["S_ts"], saved_ts_["T_ts"]
    close(saved_ts)

    saved_params_ = jldopen(joinpath(sim_output_, "model_params.jld2"))
    Sᵤ_, Sₗ_, Θᵤ_, Θₗ_ = saved_params_["Sᵤ"], saved_params_["Sₗ"],
                         saved_params_["Tᵤ"], saved_params_["Tₗ"],
    close(saved_params)

    S₀_, T₀_ = [S_ts_[:, 1, i] for i ∈ eachindex(S_ts_[1, 1, :])],
               [T_ts_[:, 1, i] for i ∈ eachindex(T_ts_[1, 1, :])]
    sal_ic_vals_string_ = string.(round.([S₀_[i][end] for i ∈ 1:num_ics]; digits = 3))

    # Salinity depth
    for i ∈ 1:num_ics
        lines!(ax[2, 1], S₀_[i], z, color = ic_colour[i, j], label = sal_ic_vals_string_[i])
    end
    # Temperature depth
    Θᵤ_array_ = fill(Θᵤ_, num_ics)
    # lines!(ax[1, 1], T₀_[1], z; label = "ΔΘ = $(abs(Θᵤ_ + Θₗ_))°C",
    #        color = T₀_[1], colorrange = Θ_range, colormap = :thermal)
    lines!(ax[1, 1], T₀_[1], z; label = "ΔΘ = $(ΔΘ_thres_vals[j])°C",
            color = ΔΘ_colour[j])
    # Salinity on Fofonoff diagram
    scatter!(ax[2, 2], [S₀_[i][end] for i ∈ 1:num_ics], Θᵤ_array_;
            color = ic_colour[:, j], markersize = 6)
    if j == 1
        axislegend(ax[2, 2], position = :lt)
    end
    Θₗ_ = [T_ts_[80, 1, i] for i ∈ 1:num_ics]
    Sₗ_ = [S_ts_[80, 1, i] for i ∈ 1:num_ics]
    Θᵤ_ = [T_ts_[81, 1, i] for i ∈ 1:num_ics]
    Sᵤ_ = [S_ts_[81, 1, i] for i ∈ 1:num_ics]

    αₗ_ = gsw_alpha.(Sₗ_, Θₗ_, p_ref)
    βₗ_ = gsw_beta.(Sₗ_, Θₗ_, p_ref)

    Δρ_thres = @. gsw_rho(Sₗ_ - (αₗ_ / βₗ_) * ΔΘ_thres_vals[j],
                          Θₗ_ - ΔΘ_thres_vals[j], p_ref) -
                  gsw_rho(Sₗ_, Θₗ_, p_ref)
    Δρ_static = @. gsw_rho(Sᵤ_, Θᵤ_, p_ref) - gsw_rho(Sₗ_, Θₗ_, p_ref)
    scatter!(ax[1, 2], fill(ΔΘ_thres_vals[j], 3), Δρ_static; color = ic_colour[:, j])
    lines!(ax[1, 2], range(ΔΘ_thres_vals[j]-0.1, ΔΘ_thres_vals[j]+0.1; length=3), Δρ_thres;
             label = "Initial ΔΘ = $(ΔΘ_thres_vals[j])",
             color = ΔΘ_colour[j])
#    if j == 3
#     axislegend(ax[1, 2]; position = :lb)
#    end
end
ic_plot

## Legend and save
#delete!(ax[1, 2])
Legend(ic_plot[3, 1], ax[1, 1], "Initial ΔΘ at the mixed\nand lower layer interface",
       orientation = :horizontal, nbanks = 3)
Legend(ic_plot[3, 2], ax[2, 1], "Initial salinity in the mixed layer",
       orientation = :horizontal, nbanks = 3)
ic_plot
colsize!(ic_plot.layout, 1, Auto(0.4))
ic_plot
#save(joinpath(PLOTDIR, "simulations", "ΔΘ_ics.png"), ic_plot)

############################################################################################
## Diffusivity time series from the simulation
############################################################################################
κ_ts_fig = Figure(resolution = (1200, 1200))
ax = [Axis(κ_ts_fig[j, i],
        xlabel = "Time (days)",
        xaxisposition = :top,
        ylabel = "Depth (metres)") for i ∈ 1:num_ics, j ∈ 1:3]
ΔΘ_vals = [0.5, 1.0, 2.0]
for (j, sim) ∈ enumerate(simulations)

    sim_output_ = joinpath(sim_datadir, sim)
    saved_ts_ = jldopen(joinpath(sim_output_, "output_timeseries.jld2"))
    S_ts_, κ_ts_ = saved_ts_["S_ts"], saved_ts_["κ_ts"]
    close(saved_ts)

    sal_ics = round.(reshape(S_ts_[end, 1, 1:end], :); digits = 4)

    for i ∈ 1:num_ics
        ax[i, j].title = "Diffusivity - initial salintiy in ML = $(sal_ics[i]),\n initial ΔΘ between layers = $(ΔΘ_vals[j])°C"
        heatmap!(ax[i, j], t, z, κ_ts_[:, :, i]', colormap = :Spectral)
    end

end
for i ∈ 4:9
    hidexdecorations!(ax[i])
end
for i ∈ [2, 3, 5, 6, 8, 9]
    hideydecorations!(ax[i])
end
tickvals = ["1e-5", "1"]
Colorbar(κ_ts_fig[:, 4];
         colormap = cgrad(:Spectral, 2, categorical = true),
         label = "Diffusivity (m²s⁻¹)", ticks = ([0.25, 0.75], tickvals))
κ_ts_fig
save(joinpath(PLOTDIR, "simulations", "ΔΘ_κ_ts.png"), κ_ts_fig)

############################################################################################
## Δρ static and cabbeling time series
############################################################################################
Δρ_s, Δρ_c = Δρ_vals["Δρ_s"], Δρ_vals["Δρ_c"]
densitydiff_TS = Figure(resolution = (1000, 700))
titles = ["Maximum static denstiy difference\nfor "*simulations[sim_num],
          "Maximum cabbeling density difference\nfor "*simulations[sim_num]]
ylabs = [L"\Delta \sigma_{0}^{s}", L"\Delta \sigma_{0}^{c}"]
ylimits = (minimum(Δρ_s)-0.002, maximum(Δρ_c)+0.001)
ax = [Axis(densitydiff_TS[1, i],
        title = titles[i],
        xlabel = "Time (days)",
        ylabel = ylabs[i],
        ylabelsize = 18) for i ∈ 1:2]
ylims!(ax[1], ylimits)
ylims!(ax[2], ylimits)
linkyaxes!(ax[1], ax[2])

line_colour = reverse(get(ColorSchemes.viridis, range(0, 1, length = num_ics)))

for i ∈ 1:num_ics
    lines!(ax[1], t, Δρ_s[:, i], color = line_colour[i])
    scatter!(ax[1], [t[1]], [Δρ_s[1, i]], color = line_colour[i],
            label = sal_ic_vals_string[i])
    #lines!(ax[1], t, Δρ_thres, color = :red)
    lines!(ax[2], t, Δρ_c[:, i], color = line_colour[i])
    scatter!(ax[2], [t[1]], [Δρ_c[1, i]], color = line_colour[i])
end

densitydiff_TS[2, :] = Legend(densitydiff_TS, ax[1],
                            "Initial salinity in the mixed layer (g/kg)";
                            orientation = :horizontal)

densitydiff_TS

############################################################################################
## Density difference threshold for time series of density difference
############################################################################################

Θₗ_ts = T_ts[80, :, 1]
Sₗ_ts = S_ts[80, :, 1]
αₗ_ts = gsw_alpha.(Sₗ_ts, Θₗ_ts, p_ref)
βₗ_ts = gsw_beta.(Sₗ_ts, Θₗ_ts, p_ref)
ΔΘ_thres_vals = [0.25, 0.5, 2.0] # set in the respective `collect_timeseires.jl` files
ΔΘ_ts = ones(length(Θₗ_ts)) .* ΔΘ_thres_vals[sim_num]
Δρ_thres_ts = @. gsw_rho(Sₗ_ts - (αₗ_ts / βₗ_ts) * ΔΘ_ts, Θₗ_ts - ΔΘ_ts, p_ref) -
                 gsw_rho(Sₗ_ts, Θₗ_ts, p_ref)

lines!(ax[1], t, Δρ_thres_ts; color = :red)
densitydiff_TS

## Threshold density diff scatter plot
fig = Figure(size = (500, 500))
ax = Axis(fig[1, 1],
          xlabel = "ΔΘ (°C)",
          ylabel = "Δρ (kgm⁻³)")
ΔΘ_thres_vals = [0.5, 1.0, 2.0]

for (j, sim) ∈ enumerate(simulations)
    sim_output_ = joinpath(SIM_DATADIR, sim)
    saved_ts_ = jldopen(joinpath(sim_output_, "output_timeseries.jld2"))
    S_ts_, T_ts_ = saved_ts_["S_ts"], saved_ts_["T_ts"]
    close(saved_ts)

    Θₗ_ = [T_ts_[80, 1, i] for i ∈ 1:num_ics]
    Sₗ_ = [S_ts_[80, 1, i] for i ∈ 1:num_ics]
    Θᵤ_ = [T_ts_[81, 1, i] for i ∈ 1:num_ics]
    Sᵤ_ = [S_ts_[81, 1, i] for i ∈ 1:num_ics]
    println(Θᵤ_)
    println(Sᵤ_)
    println(Θₗ_)
    println(Sₗ_)
    αₗ_ = gsw_alpha.(Sₗ_, Θₗ_, p_ref)
    βₗ_ = gsw_beta.(Sₗ_, Θₗ_, p_ref)

    Δρ_thres = @. gsw_rho(Sₗ_ - (αₗ_ / βₗ_) * ΔΘ_thres_vals[j], Θₗ_ - ΔΘ_thres_vals[j], p_ref) -
                  gsw_rho(Sₗ_, Θₗ_, p_ref)
    Δρ_static = @. gsw_rho(Sᵤ_, Θᵤ_, p_ref) - gsw_rho(Sₗ_, Θₗ_, p_ref)
    println(Δρ_thres)
    println(Δρ_static)
    scatter!(ax, fill(ΔΘ_thres_vals[j], 3), Δρ_static; color = ic_colour)
    lines!(ax, range(ΔΘ_thres_vals[j]-0.1, ΔΘ_thres_vals[j]+0.1; length = 3), Δρ_thres;
             label = "Δρ threshold for initial ΔΘ = $(ΔΘ_thres_vals[j])", color = :red)
end
axislegend(ax; position = :rt)
fig
