### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ 47aad13e-b752-4fac-b11b-13bcc961e057
begin
	using Pkg
	Pkg.activate("..")
	using CairoMakie, GibbsSeaWater
end

# ╔═╡ e8247e6e-3667-11ee-0293-c7e85eafeb73
md"""
# Some plots for a presentation
"""

# ╔═╡ c234a50c-848b-418b-b357-b54bd6833896
md"""
## Variables
"""

# ╔═╡ 1d1213a2-7f6e-438b-8dc6-bd43d6119369
begin
	S_star, Θ_star = 34.7, 0.5
	S_s, Θ_s = 34.55, -1.88
	N = 2000
	S_range, Θ_range = range(34.52, 34.72, length = N), range(-2, 1, length = N)
	S_grid, Θ_grid = ones(N) .* S_range', ones(N)' .* Θ_range
	ρ = gsw_rho.(S_grid, Θ_grid, 0)
	ρ_star = gsw_rho(S_star, Θ_star, 0)
	α_star = gsw_alpha(S_star, Θ_star, 0)
	β_star = gsw_beta(S_star, Θ_star, 0)
	ρ_s = gsw_rho(S_s, Θ_s, 0)
	find_Θ = findfirst(Θ_range .> -1.88)
	find_S = findfirst(ρ[find_Θ, :] .> ρ_star)
	S_iso, Θ_iso = S_range[find_S], Θ_range[find_Θ]
	gsw_rho(S_iso, Θ_iso, 0)
	nothing
end

# ╔═╡ 754e3db2-2d8c-4d14-94ea-f7bc2ba3ecd8
md"""
# Temperature and salinity plots
"""

# ╔═╡ 621b51bc-742a-487f-8e87-60a451e9c3b7
let
	fig = Figure(size = (500, 500), fontsize = 22)
	ax = Axis(fig[1, 1];
			  title = "Water masses in salinity-temperature space",
			  xlabel = "Absolute salinity (gkg⁻¹)",
			  ylabel = "Conservative temperature (°C)",
			  limits = (extrema(S_range), extrema(Θ_range)))
	contour!(ax, S_range, Θ_range, ρ'; levels = [ρ_star], labels = true, color = :black, linewidth = 0.4, labelsize = 18)
	#contour!(ax, S_range, Θ_range, ρ'; levels = [ρ_s], labels = true, color = :blue, linewidth = 0.4, labelsize = 18)
	scatter!(ax, [S_star], [Θ_star]; color = :red, label = "Deep water")
	#scatter!(ax, [S_s], [Θ_s]; color = :blue, label = "Surface water")
	scatter!(ax, [S_iso], [Θ_iso]; color = :blue, label = "Surface water")
	axislegend(ax, position = :rb)
	fig
	#save("wm_ts2.png", fig)
end

# ╔═╡ bb6d5519-9afd-4076-82bf-09c630e2b57f
let

	data = [vcat(fill(0, 5), fill(1, 5)) vcat(fill(0, 5), fill(1, 5))]
	fig = Figure(size = (500, 500), fontsize = 22)
	ax = Axis(fig[1, 1],
			  aspect = 1/3,
	          ylabel = "← z",
        	  yticksvisible = false,
        	  yticklabelsvisible = false,
        	  ygridvisible = false)
	hidexdecorations!(ax)
	hidespines!(ax)
	hm = heatmap!(ax, 1:2, 1:10, data'; colormap = [:red, :blue])
	text!(ax, 1.4, 7.5,  text = L"ρ", fontsize = 34)
	text!(ax, 1.4, 2.5,  text = L"ρ", fontsize = 34)
	fig
	#save("stable_column2.png", fig)
end

# ╔═╡ 66fc488b-34aa-4b65-8931-44dc4e72e247
let
	fig = Figure(size = (500, 500), fontsize = 22)
	ax = Axis(fig[1, 1];
			  title = "Mixing water masses in salinity-temperature space",
			  xlabel = "Absolute salinity (gkg⁻¹)",
			  ylabel = "Conservative temperature (°C)",
			  limits = (extrema(S_range), extrema(Θ_range)))
	contour!(ax, S_range, Θ_range, ρ'; levels = [ρ_star], labels = true, color = :black, linewidth = 0.4, labelsize = 18)
	#contour!(ax, S_range, Θ_range, ρ'; levels = [ρ_s], labels = true, color = :blue, linewidth = 0.4, labelsize = 18)
	scatter!(ax, [S_star], [Θ_star]; color = :red, label = "Deep water")
	#scatter!(ax, [S_s], [Θ_s]; color = :blue, label = "Surface water")
	scatter!(ax, [S_iso], [Θ_iso]; color = :blue, label = "Surface water")
	lines!(ax, [S_iso, S_star], [Θ_iso, Θ_star]; color = :purple, linestyle = :dash, label = "Mixed water")
	axislegend(ax, position = :rb)
	fig
	#save("wm_ts_mix2.png", fig)
end

# ╔═╡ 0b82cb02-7a3d-4aeb-87c1-01ea0f9f3d4f
let

	data = [vcat(fill(0, 4), 1, fill(2, 4)) vcat(fill(0, 4), 1, fill(2, 4))]
	fig = Figure(size = (500, 500), fontsize = 22)
	ax = Axis(fig[1, 1],
			  aspect = 1/3,
	          ylabel = "← z",
        	  yticksvisible = false,
        	  yticklabelsvisible = false,
        	  ygridvisible = false,
			  limits = (1, 2, 1, 10))
	hidexdecorations!(ax)
	hidespines!(ax)
	hm = heatmap!(ax, 1:2, 1:10, data'; colormap = [:red, :purple, :blue])
	text!(ax, 1.4, 7.5,  text = L"ρ", fontsize = 34)
	text!(ax, 1.4, 2.5,  text = L"ρ", fontsize = 34)
	text!(ax, 1.5, 5.55, text = L"ρ_{\mathrm{mix}}~>~ρ", fontsize = 34, align = (:center, :center))
	fig
	#save("mix_column2.png", fig)
end

# ╔═╡ c261ecbc-63d7-41fc-96ba-69f9f92c94b9
md"""
# Fofonoff diagram
"""

# ╔═╡ 274fa76d-73c5-4a82-ba78-e686e44b8119
let
	fig = Figure(size = (500, 500), fontsize = 22)
	ax = Axis(fig[1, 1];
			  title = "Fofonoff diagram",
			  xlabel = "Absolute salinity (gkg⁻¹)",
			  ylabel = "Conservative temperature (°C)",
			  limits = ((34.5, 34.72), extrema(Θ_range)))
	contour!(ax, S_range, Θ_range, ρ'; levels = [ρ_star], labels = true, color = :black, linewidth = 0.4, labelsize = 18)
	scatter!(ax, [S_star], [Θ_star]; color = :red, label = "Deep water")
	S_l = @. S_star + (α_star / β_star) * (Θ_range - Θ_star)
	lines!(ax, S_l, Θ_range; color = :red, linestyle = :dash, label = "Linearised density")
	scatter!(ax, [S_s], [Θ_s]; color = :blue, label = "Surface water")
	lines!(ax, [S_s, S_star], [Θ_s, Θ_star]; color = :purple, linestyle = :dot, label = "Mixed water")
	axislegend(ax, position = :rb)
	fig
	#save("Fof2.png", fig)
end

# ╔═╡ Cell order:
# ╟─e8247e6e-3667-11ee-0293-c7e85eafeb73
# ╟─47aad13e-b752-4fac-b11b-13bcc961e057
# ╟─c234a50c-848b-418b-b357-b54bd6833896
# ╠═1d1213a2-7f6e-438b-8dc6-bd43d6119369
# ╟─754e3db2-2d8c-4d14-94ea-f7bc2ba3ecd8
# ╟─621b51bc-742a-487f-8e87-60a451e9c3b7
# ╟─bb6d5519-9afd-4076-82bf-09c630e2b57f
# ╟─66fc488b-34aa-4b65-8931-44dc4e72e247
# ╟─0b82cb02-7a3d-4aeb-87c1-01ea0f9f3d4f
# ╟─c261ecbc-63d7-41fc-96ba-69f9f92c94b9
# ╟─274fa76d-73c5-4a82-ba78-e686e44b8119
