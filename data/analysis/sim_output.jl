using JLD2, ColorSchemes
using .VerticalProfileStability

############################################################################################
## Choose which simulation
############################################################################################
simulations = ("initial_Θ_minus1_85_0_5", "initial_Θ_18_15_20_5")
sim_num = 1 # Chnage this to look at other simulations
sim_output = joinpath(sim_datadir, simulations[sim_num])

############################################################################################
## Open and load data from chosen simulation, then extract initial conditions
############################################################################################
saved_ts = jldopen(joinpath(sim_output, "output_timeseries.jld2"))
t, S_ts, T_ts, κ_ts, Δρ_ts_0_5 = saved_ts["t"], saved_ts["S_ts"], saved_ts["T_ts"],
                                 saved_ts["κ_ts"], saved_ts["Δρ_ts"]
close(saved_ts)

## Other ΔΘ_thres values (if they have been computed)
Δρ_ts_0_25 = load(joinpath(sim_output, "ΔΘ_thres_0_25_timeseries.jld2"))["Δρ_ts"]
Δρ_ts_1 = load(joinpath(sim_output, "ΔΘ_thres_1_timeseries.jld2"))["Δρ_ts"]
Δρ_ts_2 = load(joinpath(sim_output, "ΔΘ_thres_2_timeseries.jld2"))["Δρ_ts"]

z = -500:5:-5 # This can be accessed from the saved data but much faster this way

saved_params = jldopen(joinpath(sim_output, "model_params.jld2"))
Sᵤ, Sₗ, Sᵣ, Sₘ, Θᵤ, Θₗ, num_ics = saved_params["Sᵤ"], saved_params["Sₗ"], saved_params["Sᵣ"],
                                  saved_params["Sₘ"], saved_params["Tᵤ"], saved_params["Tₗ"],
                                  saved_params["num_ics"]
if num_ics != 10; num_ics = 10; end # Can remove this if 11th sim from "initial_Θ_minus1_85_0_5" is removed.
close(saved_params)

S₀, T₀ = [S_ts[:, 1, i] for i ∈ 1:length(S_ts[1, 1, :])],
         [T_ts[:, 1, i] for i ∈ 1:length(T_ts[1, 1, :])]
sal_ic_vals_string = string.(round.([S₀[i][end] for i ∈ 1:10]; digits = 3))

############################################################################################
## Initial conditions for the simulation
############################################################################################
T = range(Θᵤ - 0.5, Θₗ + 0.5, 200)
T_grid = T' .* ones(length(T))
S = range(S₀[1][end] - 0.02, Sₘ, 200)
S_grid = S .* ones(length(S))'

ρ = gsw_rho.(S_grid, T_grid, p_ref)

lower_isopycnal = gsw_rho(Sₗ, Θₗ, p_ref)

αₗ = gsw_alpha(Sₗ, Θₗ, p_ref)
βₗ = gsw_beta(Sₗ, Θₗ, p_ref)
m = βₗ / αₗ

tang_length = 7:findfirst(S .> Sₗ)
S_tangent = S[tang_length]
tangent = @. Θₗ + m * (S_tangent - Sₗ)

Θᵤ_array = fill(Θᵤ, num_ics)
ic_colour = reverse(get(ColorSchemes.viridis, range(0, 1, length = num_ics)))

ic_plot = Figure(resolution = (1000, 1000))

xlabs = ["Salinity (g/kg)" "Temperature (∘C)"; "Salinity (g/kg)" "Temperature (∘C)"]
ylabs = ["Depth (m)" "Depth (m)"; "Temperature (∘C)" "Temperature (∘C)"]
titles = ["Salinity ICs" "Temperature ICs";
         "Fofonoff diagram - used to find the\nsalinity initial conditions in mixed alyer" "Nothing"]
ax = [Axis(ic_plot[j, i],
        xlabel = xlabs[j, i],
        xaxisposition = (j == 1 ? :top : :bottom),
        ylabel = ylabs[j, i],
        title = titles[j, i]) for i ∈ 1:2, j ∈ 1:2]
delete!(ax[2, 2])
linkyaxes!(ax[2, 1], ax[1, 1])
linkxaxes!(ax[1, 2], ax[1, 1])

# Salinity depth
for i ∈ 1:num_ics
    lines!(ax[1, 1], S₀[i], z, color = ic_colour[i], label = sal_ic_vals_string[i])
end
#Temperature depth
lines!(ax[2, 1], T₀[1], z, color = T₀[1])
ylims!(ax[2, 1], high = 0)
# Fofonoff diagram
contour!(ax[1, 2], S, T, ρ,
        color = :black,
        linewidth = 2,
        levels = [lower_isopycnal],
        label = "Isopycnal (1027.71)")
scatter!(ax[1, 2], [Sₗ], [Θₗ], color = :red, label = "Deep water mass")
lines!(ax[1, 2], S_tangent, tangent, color = :red, label = "Linearised density about\ndeep water mass")
#lines!(ax[1, 2], Sₗ:0.05:Sₘ, 0.5 * ones(length(Sₗ:0.05:Sₘ)), color = :orange, linestyle = :dot, label = "Salinity in lower layer")
scatter!(ax[1, 2], [S₀[i][end] for i ∈ 1:10], Θᵤ_array;
         color = ic_colour, markersize = 4)
axislegend(ax[1, 2], position = :lt)
ic_plot[2, 2] = Legend(ic_plot, ax[1, 1], "Initial salinity in the mixed layer")
ic_plot
#save(joinpath(plotdir, "simulations", simulations[sim_num], "ics.png"), ic_plot)

############################################################################################
## Diffusivity time series from the simulation
############################################################################################
κ_ts_fig = Figure(resolution = (1400, 800))
sal_ics = round.([reshape(S_ts[end, 1, 1:5], :) reshape(S_ts[end, 1, 6:10], :)]';
                 digits = 3)
ax = [Axis(κ_ts_fig[i, j],
           title = "Diffusivity - salintiy = $(sal_ics[i, j])",
           xlabel = "Time (days)",
           xaxisposition = :top,
           ylabel = "Depth (metres)") for i ∈ 1:2, j ∈ 1:5]
k = 1
for i ∈ 1:2, j ∈ 1:5
    heatmap!(ax[i, j], t, z, κ_ts[:, :, k]', colormap = :Spectral)
    k += 1
end
tickvals = ["1e-5", "1"]
Colorbar(κ_ts_fig[:, 6];
         colormap = cgrad(:Spectral, 2, categorical = true),
         label = "Diffusivity (m²s⁻¹)", ticks = ([0.25, 0.75], tickvals))
κ_ts_fig
#save(joinpath(plotdir, "simulations", simulations[sim_num], "κ_ts.png"), κ_ts_fig)

############################################################################################
## Δρ static and cabbeling time series
############################################################################################
density_diffs = [Δρ_ts_0_25, Δρ_ts_0_5, Δρ_ts_1, Δρ_ts_2]
ΔΘ_thres_val = [0.25, 0.5, 1.0, 2.0]
ΔΘ_thres_idx = 4
Δρ_s, Δρ_c = density_diffs[ΔΘ_thres_idx]["Δρ_s"], density_diffs[ΔΘ_thres_idx]["Δρ_c"]
densitydiff_TS = Figure(resolution = (1000, 700))
titles = ["Maximum static denstiy difference\nfor ΔΘ_thres = $(ΔΘ_thres_val[ΔΘ_thres_idx])ᵒC",
          "Maximum cabbeling density difference\nfor ΔΘ_thres = $(ΔΘ_thres_val[ΔΘ_thres_idx])ᵒC"]
ylabs = [L"\Delta \sigma_{0}^{s}", L"\Delta \sigma_{0}^{c}"]
ylimits = (minimum(Δρ_s)-0.002, maximum(Δρ_c))
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
#save(joinpath(plotdir, "simulations", simulations[sim_num], "Δρ_ts_0_225.png"), densitydiff_TS)

############################################################################################
## Density difference threshold for time series of density difference
############################################################################################
Θₗ_ts = T_ts[80, :, 10]
Sₗ_ts = S_ts[80, :, 10]
αₗ_ts = gsw_alpha.(Sₗ_ts, Θₗ_ts, p_ref)
βₗ_ts = gsw_beta.(Sₗ_ts, Θₗ_ts, p_ref)
ΔΘ_ts = ones(length(Θₗ_ts)) .* ΔΘ_thres_val[ΔΘ_thres_idx]
Δρ_thres_ts = @. gsw_rho(Sₗ_ts - (αₗ_ts / βₗ_ts) * ΔΘ_ts, Θₗ_ts - ΔΘ_ts, p_ref) -
                 gsw_rho(Sₗ_ts, Θₗ_ts, p_ref)

lines!(ax[1], t, Δρ_thres_ts; color = :red)
densitydiff_TS
