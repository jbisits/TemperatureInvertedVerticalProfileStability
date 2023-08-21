### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ 07a40b83-ddc8-4c90-94b0-a1c0780e8fc8
begin
	using Pkg
	Pkg.activate("..")
end

# ╔═╡ 76f0d455-f0ea-4465-a86c-548cdcdff11c
using PlutoUI

# ╔═╡ cd4566fb-c16c-4c4c-b308-5203dae44acd
begin
	ECCO_plots_path = joinpath(@__DIR__, "..", "plots", "ECCO")
	ECCO_plots = "../plots/ECCO/" .* readdir(ECCO_plots_path)
	nothing
end

# ╔═╡ 32ef1b82-a083-11ed-0143-b7d663e4f60e
md"""
# ECCO data output

This notebook contains plots of analysis from the ECCO data.
We are looking at the maximum density difference of a profile against the temperature of the lower level.
The static density difference is calculated for all pairs of levels that have a temperature difference greater than some chosen threshold, and the maximum is found to give an upper and lower level.
The cabbeling density difference is computed between the upper and lower level.
The plots below are of the maximum static (left panel) and cabbeling (right panel) density difference of a profile against the temperature of the *lower level*.

Varying the temperature threshold allows us to see what happens when there are larger temperature differences between the levels increasing the cabbeling instability.

For ``\Delta\Theta`` I have looked at so far there are three plots.
A full plot, then two versions of this same plot zoomed in on the y-axis (``\Delta\rho``) to see what is happening closer to the instability.

**Note:** this notebook does not contain any code (yet) I just thought it would be a good way to show plots.

## ``\Delta\Theta = 0.5^{\circ}C``

### Full plot
$(LocalResource(ECCO_plots[1]))
### Zoom 1
$(LocalResource(ECCO_plots[2]))
### Zoom 2
$(LocalResource(ECCO_plots[3]))

## ``\Delta\Theta = 1.0^{\circ}C``

### Full plot
$(LocalResource(ECCO_plots[4]))
### Zoom 1
$(LocalResource(ECCO_plots[5]))
### Zoom 2
$(LocalResource(ECCO_plots[6]))

## ``\Delta\Theta = 2.0^{\circ}C``

### Full plot
$(LocalResource(ECCO_plots[7]))
### Zoom 1
$(LocalResource(ECCO_plots[8]))
### Zoom 2
$(LocalResource(ECCO_plots[9]))
"""

# ╔═╡ 0710fc63-729d-469e-a63e-f71032f182d4
TableOfContents(title = "ECCO plots")

# ╔═╡ Cell order:
# ╟─07a40b83-ddc8-4c90-94b0-a1c0780e8fc8
# ╟─76f0d455-f0ea-4465-a86c-548cdcdff11c
# ╟─cd4566fb-c16c-4c4c-b308-5203dae44acd
# ╟─32ef1b82-a083-11ed-0143-b7d663e4f60e
# ╟─0710fc63-729d-469e-a63e-f71032f182d4
