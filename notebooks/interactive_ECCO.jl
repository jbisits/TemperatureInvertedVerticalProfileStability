### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ a4d1faa4-f76d-49fe-99cf-f21ea0d1ee29
begin
	using Pkg
	Pkg.activate("..")
	using Rasters, NCDatasets, OceanRasterConversions, PlutoUI, CairoMakie, Glob, Dates, GibbsSeaWater
end

# ╔═╡ 70810fe8-3fc1-11ee-3fca-8bbee662a36b
md"""

# Reactivity

[Pluto.jl](https://github.com/fonsp/Pluto.jl) is a Julia package for *simple, reactive notebooks*.

Similar to Jupyter, the notebooks are structured around cells that can contain markdown or code.

The difference is that when (global) variables are changed in a cell, **Pluto.jl will re-execute all the cells that contain this variable.**
"""

# ╔═╡ 49e79b9f-bfcc-478e-90ff-6b1ee9ed5344
N = -10

# ╔═╡ 1101a1bb-5640-4459-86d6-3b8ff5b96ce6
N^2

# ╔═╡ ec5d0761-7cdf-47bf-b5f2-456a8ac6a7ea
let
	x = range(-5, 5, length = 250)
	y = @. x^2 - N
	sign_string = sign(N) == 1 ? "-" : "+" 
	fig, ax = lines(x, y; label = "y = x² $(sign_string) $(abs(N))")
	ax.xlabel = "x"
	ax.ylabel = "y"
	axislegend(ax, position = :rb)
	fig
end

# ╔═╡ cc31a3d7-6279-4986-8c49-b33b22ebff12
begin
	ECCO_2007_folder = joinpath(@__DIR__, "../data/observations/ECCO_daily_mean_TS")
	ECCO_files = glob("*.nc", ECCO_2007_folder)
	days = Date(2007, 01, 01):Day(1):Date(2007, 12, 31)
	day_slider = @bind day PlutoUI.Slider(days, show_value = true)
	r = Raster(ECCO_files[1], name = :SALT)
	x, y, z = lookup(r, X), lookup(r, Y), lookup(r, Z)
	long_choose = @bind long TextField(default = "0")
	#long_choose = @bind long Select(x)
	lat_choose = @bind lat TextField(default = "-60")
	#lat_choose = @bind lat Select(y)
	depth_choose = @bind depth Select(z)
	nothing
end

# ╔═╡ 4af9712a-de33-47ba-80c0-a1436dacf9a5
md"""
# Interactivity

The interactive side of Pluto.jl comes from [PlutoUI.jl](https://github.com/JuliaPluto/PlutoUI.jl) which makes it easy to use `html` commands in Pluto.jl to make these notebooks *reactive*.
"""

# ╔═╡ a8e22201-c486-4bc9-932e-4c24ff39ee27
md"""
# Exploring ECCOv4r4 data interactively

Looking at the daily salinity and temperature data from 2007.
Using Julia and [Pluto.jl](https://github.com/fonsp/Pluto.jl) we can view data in an easy to setup *interactive* manner.

## Temperature and salinity plots at depth

Date = $day_slider

Depth = $depth_choose m
"""

# ╔═╡ b6b470c7-51d6-41b5-a7e0-e0d8e3a51d60
begin	
	day_idx = findfirst(day .== days)
	depth_idx = findfirst(z .== depth)
	rs = RasterStack(ECCO_files[day_idx], name = (:SALT, :THETA))
	vars = keys(rs)
	fig = Figure(size = (800, 800))
	ax = [Axis(fig[i, 1]; xlabel = "Longitude", ylabel = "Latitude") for i ∈ eachindex(vars)]
	for (i, var) ∈ enumerate(vars)
		colormap = var == :SALT ? :haline : :thermal
		hm = heatmap!(ax[i], rs[var][Z(depth_idx)]; colormap)
		ax[i].title = string(var) * ", z = $depth m, $day."
		unit = var == :SALT ? "psu" : "°C"
		Colorbar(fig[i, 2], hm, label = unit)
	end
	hidexdecorations!(ax[1], grid = false)
	fig
end

# ╔═╡ 38561f85-9421-4602-8b8d-1af8d9048964
md"""

## Temperature and salinity profiles

Longitude = $long_choose

Latitude = $lat_choose
"""

# ╔═╡ b3d09471-0c8b-436a-9d1b-87383dc33240
let
	day_idx = findfirst(day .== days)
	long = parse(Float32, long)
	lat = parse(Float32, lat)
	rs = RasterStack(ECCO_files[day_idx], name = (:SALT, :THETA))
	S = vec(rs[:SALT][X(Near(long)), Y(Near(lat))].data)
	T = vec(rs[:THETA][X(Near(long)), Y(Near(lat))].data)
	z = lookup(rs, Z)
	fig, ax, plt = scatter(S, T, color = z)
	Colorbar(fig[1, 2], plt)
	ax.xlabel = "S (psu)"
	ax.ylabel = "θ (°C)"
	ax.title = "T-S profile at $(abs(long))°$(ifelse(long ≥ 0, "E", "W")),  $(abs(lat)) ° $(ifelse(lat ≤ 0, "S", "N")), $day"

	Θᵤ =0.073116146
	upper_idx = findfirst(T .≥ Θᵤ)
	Θₗ = 1.1258758
	lower_idx = isnothing(findfirst(skipmissing(T) .≥ Θₗ)) ? findmax(skipmissing(T))[2] : findfirst(skipmissing(T) .≥ Θₗ)
	scatter!(ax, [S[upper_idx], S[lower_idx]], [T[upper_idx], T[lower_idx]], markersize = 5)

	zᵤ = z[upper_idx]
	zₗ = z[lower_idx]
	ẑ = -(zᵤ + zₗ) / 2

	ρᵤ = gsw_rho(S[upper_idx], T[upper_idx], ẑ)
	ρₗ = gsw_rho(S[lower_idx], T[lower_idx], ẑ)
	ρᵤ - ρₗ

	S_linear_range = range(minimum(skipmissing(S)), S[lower_idx], length = 50)
	αₗ = gsw_alpha(S[lower_idx], T[lower_idx], ẑ)
	βₗ = gsw_beta(S[lower_idx], T[lower_idx], ẑ)
	Θ_linear = @. T[lower_idx] + (βₗ / αₗ) * (S_linear_range - S[lower_idx])
	lines!(ax, S_linear_range, Θ_linear)

	S_range = range(extrema(skipmissing(S))..., length = 1000)
	T_range = range(extrema(skipmissing(T))..., length = 1000)

	S_grid = S_range .* ones(length(S_range))'
	T_grid = T_range' .* ones(length(T_range))

	ρ = gsw_rho.(S_grid, T_grid, ẑ)

	contour!(ax, S_range, T_range, ρ, levels = [ρₗ])
	
	fig
	#save("Tinv_profile_Weddel.png", fig)
end

# ╔═╡ 256fc23b-d784-48da-adb4-1c4499359010
#####################################################################################
## Possible S-T profile to demonstrate
#####################################################################################
let
	long = -43
	lat = -70
	rs = RasterStack(ECCO_files[100], name = (:SALT, :THETA))
	S = vec(rs[:SALT][X(Near(long)), Y(Near(lat))].data)
	T = vec(rs[:THETA][X(Near(long)), Y(Near(lat))].data)
	
	## Linear criteria
	T_max, T_find = findmax(skipmissing(T))
	α = gsw_alpha(S[T_find], T[T_find], 0)
	β = gsw_beta(S[T_find], T[T_find], 0)
	S_linear_range = range(34.49, S[T_find], length = 20)
	T_linear = @. T[T_find] + (β / α) * (S_linear_range - S[T_find])
	
	## Density
	S_range = range(extrema(skipmissing(S))..., length = 1000)
	T_range = range(extrema(skipmissing(T))..., length = 1000)
	S_grid = S_range .* ones(length(S_range))'
	T_grid = T_range' .* ones(length(T_range))
	p_ref = 0.0 #reference pressure
	ρ = gsw_rho.(S_grid, T_grid, p_ref)
	isopycnal = gsw_rho(S[T_find], T[T_find], p_ref)
	
	fig, ax = lines(S, T, color = (:gray, 0.5), label = "ECCO profiles")
	xlims!(ax, 34.15, nothing)
	ax.xlabel = "S (psu)"
	ax.ylabel = "θ (°C)"
	ax.title = "Monthly S-T profiles at $(abs(long))°$(ifelse(long ≥ 0, "E", "W")),  $(abs(lat))°$(ifelse(lat ≤ 0, "S", "N"))"

	for d ∈ days[32:30:end] #days[214:243]
		i = findfirst(d .== days)
	    rs = RasterStack(ECCO_files[i], name = (:SALT, :THETA))
	    S_ = vec(rs[:SALT][X(Near(long)), Y(Near(lat))].data)
	    T_ = vec(rs[:THETA][X(Near(long)), Y(Near(lat))].data)
	    lines!(ax, S_, T_, color = (:gray, 0.5))
	end
	lines!(ax, S_linear_range, T_linear, label = "Cabbeling instability threshold", color = :red, linestyle = :dash)
	contour!(ax, S_range, T_range, ρ, levels = [isopycnal], label = "Static instability threshold", color = :red)
	S_synthetic = zeros(length(S[1:T_find]))
	S_synthetic[1:5] .= S[1:5]
	S_synthetic[6:17] .= S[6:17] .+ 0.06 * abs.(range(-3.1, -0.8, length = 12))
	S_synthetic[16:end] .= range(S_synthetic[15], S[T_find], length = 7)
	T_synthetic = T[1:T_find]
	T_synthetic[16:end] .= range(T_synthetic[15], T[T_find], length = 7)
	lines!(ax, S_synthetic, T_synthetic, color = :magenta, label = "Synthetic unstable to cabbeling profile")
	axislegend(ax, position = :lt)
	fig
end

# ╔═╡ cab5981e-e729-4e6b-aea3-894fd128d6f2
TableOfContents(title = "Interactivity in Julia")

# ╔═╡ Cell order:
# ╟─a4d1faa4-f76d-49fe-99cf-f21ea0d1ee29
# ╟─70810fe8-3fc1-11ee-3fca-8bbee662a36b
# ╠═49e79b9f-bfcc-478e-90ff-6b1ee9ed5344
# ╠═1101a1bb-5640-4459-86d6-3b8ff5b96ce6
# ╟─ec5d0761-7cdf-47bf-b5f2-456a8ac6a7ea
# ╠═cc31a3d7-6279-4986-8c49-b33b22ebff12
# ╟─4af9712a-de33-47ba-80c0-a1436dacf9a5
# ╟─a8e22201-c486-4bc9-932e-4c24ff39ee27
# ╠═b6b470c7-51d6-41b5-a7e0-e0d8e3a51d60
# ╟─38561f85-9421-4602-8b8d-1af8d9048964
# ╠═b3d09471-0c8b-436a-9d1b-87383dc33240
# ╠═256fc23b-d784-48da-adb4-1c4499359010
# ╟─cab5981e-e729-4e6b-aea3-894fd128d6f2
