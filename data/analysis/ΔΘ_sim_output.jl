using JLD2, ColorSchemes
using .VerticalProfileStability

############################################################################################
## Choose which simulation
############################################################################################
simulations = ("initial_ΔΘ_0.5", "initial_ΔΘ_1.0", "initial_ΔΘ_2.0") .* "_tgrad"
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
ΔΘ_thres_vals = [0.5, 1.0, 2.0]
############################################################################################
## Initial conditions for the simulation
############################################################################################

# Axis and constant setups

##
Θ = range(-2, 3, 1000)
Θ_grid = Θ' .* ones(length(Θ))
S = range(34.54, 34.75, 1000)
S_grid = S .* ones(length(S))'

p_ref = 0.0 #reference pressure
ρ = gsw_rho.(S_grid, Θ_grid, p_ref)

Sₗ, Θₗ = 34.7,  0.5
lower_isopycnal = gsw_rho(Sₗ, Θₗ, p_ref)
αₗ = gsw_alpha(Sₗ, Θₗ, p_ref)
βₗ = gsw_beta(Sₗ, Θₗ, p_ref)
m = βₗ / αₗ
find_iso = findall(ρ .≈ lower_isopycnal)
iso_S = S_grid[find_iso] # salinity values for the isopycnal
iso_Θ = Θ_grid[find_iso] # Θ values for the isopycnal
S_linear = range(34.514, S[end]; length = length(iso_S))
Θ_linear = @. Θₗ + m * (S_linear - Sₗ)

ic_colour_stability = get(ColorSchemes.dense, range(0.25, 1, length = 3))
ic_colour = get(ColorSchemes.haline, range(0, 0.95, length = 3))
ic_colour = reverse(get(ColorSchemes.haline, range(0, 0.95, length = num_ics * 3)))
ic_colour = reshape(ic_colour, 2, 3)
#ΔΘ_colour = [:blue, :orange, :green]
ΔΘ_colour = get(ColorSchemes.thermal, range(0, 0.8, length = 4))
Θ_range = (-2, 2)
expt_ls = (:dot, :dash, :dashdot)

##
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
# contour!(ax[2, 2], S, Θ, ρ;
#         label = L"Isopycnal through $(S^{*},~\Theta^{*})$",
#         color = density_grad[2],
#         linewidth = 2,
#         levels = [lower_isopycnal])
lines!(ax[2, 2], iso_S, iso_Θ;
        label = L"Isopycnal through $(S_{l},~\Theta_{l})$",
        color = :black,
        linewidth = 2)
lines!(ax[2, 2], S_linear, Θ_linear; color = :grey, linewidth = 2,
       label = L"Linearised density at $(S_{l},~\Theta_{l})$")
scatter!(ax[2, 2], [Sₗ], [Θₗ], color = :red,
        label = L"(S_{l},~\Theta_{l})")

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
        lines!(ax[2, 1], S₀_[i], z; color = ic_colour[3-i, j], linestyle = expt_ls[j],
               label = sal_ic_vals_string_[i], linewidth = 2)
    end
    # Temperature depth
    Θᵤ_array_ = fill(Θᵤ_, num_ics)
    # lines!(ax[1, 1], T₀_[1], z; label = "ΔΘ = $(abs(Θᵤ_ + Θₗ_))°C",
    #        colormap = :thermal, color = T₀_[1],
    #        linestyle = expt_ls[j], linewidth = 3)
    lines!(ax[1, 1], T₀_[1], z; label = "ΔΘ = $(ΔΘ_thres_vals[j])°C",
            linestyle = expt_ls[j], linewidth = 2, color = :black)
    # Salinity on Fofonoff diagram
    scatter!(ax[2, 2], [S₀_[i][end] for i ∈ 1:num_ics], Θᵤ_array_;
            color = reverse(ic_colour[:, j]), markersize = 6)
    # if j == 1
    #     axislegend(ax[2, 2], position = :lt)
    # end
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
    # lines!(ax[1, 2], range(ΔΘ_thres_vals[j]-0.1, ΔΘ_thres_vals[j]+0.1; length=2), Δρ_thres;
    #       label = "Δρ threshold for initial ΔΘ = $(ΔΘ_thres_vals[j])",
    #       color = ΔΘ_colour[j])
    hlines!(ax[1, 2], Δρ_thres; linestyle = expt_ls[j], color = :black,
            label = "initial ΔΘ = $(ΔΘ_thres_vals[j])")
    scatter!(ax[1, 2], fill(ΔΘ_thres_vals[j], 2), Δρ_static;
             color = reverse(ic_colour[:, j]))
   if j == 3
    #axislegend(ax[1, 1]; position = :lb, orientation = :horizontal, nbanks = 3)
    hlines!(ax[1, 2], 0;
           color = ic_colour_stability[3]#=, label = "Static instability threshold"=#)
    #axislegend(ax[1, 2]; position = :lb, orientation = :horizontal, nbanks = 4)
    axislegend(ax[1, 2], ax[1, 2], "Δρ threshold for"; position = :lb, orientation = :horizontal, nbanks = 4)
    #axislegend(ax[2, 1]; position = :lb, orientation = :horizontal, nbanks = 2)
    axislegend(ax[2, 2]; position = :lt)
   end
end
colsize!(ic_plot.layout, 1, Auto(0.6))
ic_plot
#save(joinpath(PLOTDIR, "simulations", "model_ics.png"), ic_plot)

## Legend and save
# delete!(ax[1, 2])
# Legend(ic_plot[3, 1], ax[1, 1], "Initial ΔΘ at the mixed\nand lower layer interface",
#        orientation = :horizontal, nbanks = 3)
# Legend(ic_plot[3, 2], ax[2, 1], "Initial salinity in the mixed layer",
#        orientation = :horizontal, nbanks = 3)
# ic_plot

############################################################################################
## Diffusivity time series from the simulation
############################################################################################
κ_ts_fig = Figure(resolution = (850, 1200))
ax = [Axis(κ_ts_fig[j, i],
        xlabel = "Time (days)",
        xaxisposition = :top,
        ylabel = "Depth (metres)") for i ∈ 1:num_ics, j ∈ 1:3]
ΔΘ_vals = [0.5, 1.0, 2.0]
for (j, sim) ∈ enumerate(simulations)

    sim_output_ = joinpath(SIM_DATADIR, sim)
    saved_ts_ = jldopen(joinpath(sim_output_, "output_timeseries.jld2"))
    S_ts_, κ_ts_ = saved_ts_["S_ts"], saved_ts_["κ_ts"]
    close(saved_ts)

    sal_ics = round.(reshape(S_ts_[end, 1, 1:end], :); digits = 4)

    for i ∈ 1:num_ics
        ax[i, j].title = "Diffusivity - initial salintiy in ML = $(sal_ics[i]),\n initial ΔΘ between layers = $(ΔΘ_vals[j])°C"
        heatmap!(ax[i, j], t, z, κ_ts_[:, :, i]'; colormap = :Spectral,
                 colorrange = (1e-5, 1))
    end

end
tickvals = ["1e-5", "1"]
Colorbar(κ_ts_fig[:, 3];
         colormap = cgrad(:Spectral, 2, categorical = true),
         label = "Diffusivity (m²s⁻¹)", ticks = ([0.25, 0.75], tickvals))
κ_ts_fig
#save(joinpath(PLOTDIR, "simulations", "model_diff.png"), κ_ts_fig)

## Average diffusivity over depth
fig = Figure(size = (500, 500))
ax = Axis(fig[1, 1];
          title = "Average diffusivity over below mixed layer",
          xlabel = "time",
          ylabel = "Average diffusivity below mixed layer")

for (j, sim) ∈ enumerate(simulations)

    sim_output_ = joinpath(SIM_DATADIR, sim)
    saved_ts_ = jldopen(joinpath(sim_output_, "output_timeseries.jld2"))
    S_ts_, κ_ts_ = saved_ts_["S_ts"], saved_ts_["κ_ts"]
    close(saved_ts)

    sal_ics = round.(reshape(S_ts_[end, 1, 1:end], :); digits = 4)
    m1 = vec(mean(κ_ts_[1:82, :, 1]; dims = 1))
    m2 = vec(mean(κ_ts_[1:82, :, 2]; dims = 1))
    lines!(ax, t, m1; color = ic_colour[1], label = "stable")
    lines!(ax, t, m2; color = ic_colour[2], label = "cabbeling unstable")
end
fig
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
