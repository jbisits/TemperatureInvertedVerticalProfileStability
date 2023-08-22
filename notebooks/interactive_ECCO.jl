### A Pluto.jl notebook ###
# v0.19.27

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
	using Rasters, NCDatasets, OceanRasterConversions, PlutoUI, CairoMakie, Glob, Dates
end

# ╔═╡ 70810fe8-3fc1-11ee-3fca-8bbee662a36b
md"""

# Reactivity

[Pluto.jl](https://github.com/fonsp/Pluto.jl) is a Julia package for *simple, reactive notebooks*.

Similar to Jupyter, the notebooks are structured around cells that can contain markdown or code.

The difference is that when (global) variables are changed in a cell, **Pluto.jl will re-execute all the cells that contain this variable.**
"""

# ╔═╡ 49e79b9f-bfcc-478e-90ff-6b1ee9ed5344
N = 20

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
	fig, ax = lines(S, T)
	ax.xlabel = "S (psu)"
	ax.ylabel = "θ (°C)"
	ax.title = "T-S profile at $(abs(long)) °$(ifelse(long ≤ 0, "E", "W")),  $(abs(lat)) ° $(ifelse(lat ≤ 0, "S", "N")), $day"
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
# ╟─cc31a3d7-6279-4986-8c49-b33b22ebff12
# ╟─4af9712a-de33-47ba-80c0-a1436dacf9a5
# ╟─a8e22201-c486-4bc9-932e-4c24ff39ee27
# ╟─b6b470c7-51d6-41b5-a7e0-e0d8e3a51d60
# ╟─38561f85-9421-4602-8b8d-1af8d9048964
# ╟─b3d09471-0c8b-436a-9d1b-87383dc33240
# ╟─cab5981e-e729-4e6b-aea3-894fd128d6f2
