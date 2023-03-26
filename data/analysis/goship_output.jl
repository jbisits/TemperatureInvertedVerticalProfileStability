using .VerticalProfileStability
using StatsBase, LinearAlgebra, ColorSchemes, GeoMakie

const GOSHIP_DATA = joinpath(@__DIR__, "goship.jld2")
const GOSHIP_JOINED = joinpath(@__DIR__, "goship_joined.jld2")

gd = jldopen(GOSHIP_DATA)

ΔΘ_vals = (0.5, 1.0, 2.0, 3.0)
oceans = ("atlantic", "indian", "pacific", "southern")

## Location of temperature inverted profiles (or just profiles).
# Location of inverted profiles is nice because easy to show polar regions are key.
fig = Figure(size = (500, 500))
ax = GeoAxis(fig[1, 1];
             title = "CTD ship data",
             xlabel = "Longitude",
             ylabel = "Latitude",
             coastlines = true)

for key ∈ keys(gd["1.0"])

    Θᵤ = collect(skipmissing(gd["1.0"][key]["Θᵤ"]))
    Θₗ = collect(skipmissing(gd["1.0"][key]["Θₗ"]))
    lons = gd["1.0"][key]["lons"]
    lats = gd["1.0"][key]["lats"]
    sc = scatter!(ax, lons, lats; color = :grey, markersize = 5)
    find_inv = findall(Θᵤ .< Θₗ)
    lons_inv = lons[find_inv]
    lats_inv = lats[find_inv]
    scinv = scatter!(ax, lons_inv, lats_inv; color = (:red, 0.5), markersize = 5)
    Legend(fig[2, 1], [sc, scinv], ["CTD measurement", "Temperature inverted"],
           orientation = :horizontal)

end
#Legend(fig[2, 1], ax, orientation = :horizontal)
fig
save(joinpath(PLOTDIR, "GOSHIP", "CTD_location.png"), fig)

##
fig2 = Figure()
ax2 = Axis(fig2[1, 1];
          xlabel = "Θₗ",
          ylabel = "Δρ_static")
for ΔΘ_key ∈ keys(gd)
    for key ∈ keys(gd["1.0"])
        Δρˢ = collect(skipmissing(gd[ΔΘ_key][key]["Δρˢ"]))
        Θᵤ = collect(skipmissing(gd[ΔΘ_key][key]["Θᵤ"]))
        Θₗ = collect(skipmissing(gd[ΔΘ_key][key]["Θₗ"]))
        find_inv = Θᵤ .< Θₗ
        scatter!(ax2, Θₗ[find_inv], Δρˢ[find_inv]; label = key)
    end
end
Legend(fig2[2, 1], ax, orientation = :horizontal)
fig2
ylims!(ax2, -0.1, 0.01)
xlims!(ax2, -1.8, 15)
fig2
##
close(gd)

## Joined data
gdj = jldopen(GOSHIP_JOINED)
gdj["0.5"]

ΔΘ_thres = (0.5, 1.0, 2.0, 3.0)
ΔΘ_colours = get(ColorSchemes.thermal, range(0, 0.8; length = 4))

## Figure
Θₗ_lims = (-1.88, 10)
Θₗ_range = range(Θₗ_lims...; length = 100)
Δρ_lims = (-0.2, 0.02)
bin_width = 0.001
full_fig = Figure(size = (800, 1000))
# scatter
splot = full_fig[1, 1] = GridLayout()
ax_splot = Axis(splot[1, 1];
                xlabel = "Θ (ᵒC)",
                xaxisposition = :top,
                title = "(a) Maximum static density difference between two vertically\nspaced levels of a profile against temperature of lower level",
                ylabel = "Δρ (kgm⁻³)")
data_count = Vector{Int64}(undef, length(ΔΘ_thres))
for (i, key) ∈ enumerate(keys(gdj))

    Θᵤ = collect(skipmissing(gdj[key]["Θᵤ"]))
    Θₗ = collect(skipmissing(gdj[key]["Θₗ"]))
    find_inv = Θᵤ .≤ Θₗ
    Θₗ = Θₗ[find_inv]
    Δρˢ = collect(skipmissing(gdj[key]["Δρˢ"]))[find_inv]
    data_count[i] = length(Δρˢ)

    scatter!(ax_splot, Θₗ, Δρˢ; color = ΔΘ_colours[i], markersize = 5)

    Sₗ_mean = mean(collect(skipmissing(gdj[key]["Sₗ"]))[find_inv])
    pₘ = 0.5 .* (collect(skipmissing(gdj[key]["pₗ"]))[find_inv] .+
                    collect(skipmissing(gdj[key]["pᵤ"]))[find_inv])
    pₘ_mean = mean(pₘ)
    αₗ = gsw_alpha.(Sₗ_mean, Θₗ_range, pₘ_mean)
    βₗ = gsw_beta.(Sₗ_mean, Θₗ_range, pₘ_mean)
    Δρ_thres = @. gsw_rho.(Sₗ_mean - (αₗ / βₗ) * ΔΘ_thres[i], Θₗ_range - ΔΘ_thres[i],
                            pₘ_mean) - gsw_rho.(Sₗ_mean, Θₗ_range, pₘ_mean)
    lines!(ax_splot, Θₗ_range, Δρ_thres; color = ΔΘ_colours[i], linewidth = 2,
            label = "ΔΘ = $(ΔΘ_thres[i])ᵒC")
end
full_fig
xlims!(ax_splot, Θₗ_lims)
ylims!(ax_splot, Δρ_lims)
Legend(splot[1, 2], ax_splot, "Δρ threshold for")
full_fig
data_count

# pdf
pdf_plots = full_fig[2:3, 1] = GridLayout()
ax_pdf = [Axis(pdf_plots[i, j];
        xlabel = "Δρ (kgm⁻³)")
      for i ∈ 1:2, j ∈ 1:2]
less_thres = Vector{Float64}(undef, 4)
over_thres = Vector{Float64}(undef, 4)
data_count2 = Vector{Int64}(undef, length(ΔΘ_thres))
letter_labels = ["(b)", "(c)", "(d)", "(e)"]

for (i, key) ∈ enumerate(keys(gdj))
    Θₗ = collect(skipmissing(gdj[key]["Θₗ"]))
    Θᵤ = collect(skipmissing(gdj[key]["Θᵤ"]))
    find_inverted = Θᵤ .< Θₗ
    Δρˢ = collect(skipmissing(gdj[key]["Δρˢ"]))[find_inverted]
    data_count2[i] = length(Δρˢ)
    Sₗ_mean = mean(collect(skipmissing(gdj[key]["Sₗ"]))[find_inverted])
    pₘ = 0.5 .* (collect(skipmissing(gdj[key]["pₗ"]))[find_inverted] .+
                collect(skipmissing(gdj[key]["pᵤ"]))[find_inverted])
    pₘ_mean = mean(pₘ)
    αₗ = gsw_alpha.(Sₗ_mean, Θₗ_range, pₘ_mean)
    βₗ = gsw_beta.(Sₗ_mean, Θₗ_range, pₘ_mean)
    Δρ_thres = @. gsw_rho.(Sₗ_mean - (αₗ / βₗ) * ΔΘ_thres[i], Θₗ_range - ΔΘ_thres[i],
                        pₘ_mean) - gsw_rho.(Sₗ_mean, Θₗ_range, pₘ_mean)
    Δρ_thres_mean = mean(Δρ_thres)
    hist_edges = minimum(Δρˢ):bin_width:maximum(Δρˢ)

    hist!(ax_pdf[i], Δρˢ; bins = hist_edges, normalization = :pdf,
          color = (ΔΘ_colours[i], 0.5))
    vlines!(ax_pdf[i], Δρ_thres_mean; color = ΔΘ_colours[i], linewidth = 2,
            label = "Δρ threshold for ΔΘ")
    vlines!(ax_pdf[i], 0; color = :black, linestyle = :dash,
            label = "Static stability threshold")
    ax_pdf[i].title = letter_labels[i] * " PDF for ΔΘ = $(ΔΘ_thres[i])°C"
    hist_fit = fit(Histogram, Δρˢ, hist_edges)
    hist_fit = normalize(hist_fit; mode = :pdf)

    find_thres = findall(hist_edges .≤ Δρ_thres_mean)

    # to average threshold
    less_thres[i] = sum(hist_fit.weights[find_thres] .* bin_width)
    # after average threshold
    over_thres[i] = 1 - less_thres[i]

end
full_fig
pdf_lims = (-0.2, 0.01)
for i ∈ 1:4
    xlims!(ax_pdf[i], pdf_lims)
end
less_thres
over_thres
rowsize!(full_fig.layout, 1, Auto(1.15))
full_fig
data_count == data_count2
data_count
#save(joinpath(PLOTDIR, "GOSHIP", "goship_sc_pdf.png"), full_fig)
##
close(gdj)
