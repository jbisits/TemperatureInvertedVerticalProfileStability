# Vertical profile stability

Look at the stability of vertical profiles in the ocean to assess if cabbeling shapes the nature of profiles in the global ocean.

To reproduce the results in the repository first clone the repo and change to the folder where it is saved.
The project then needs to be activated and instantiated (build all required dependencies) to be used,

```julia
julia> ] # to enter the package manager
(@v1.8)> activate .
(VerticalProfileStability)> instantiate
```

After this all code in the repository should be able to be executed and results viewed (provided the correct data has been downloaded).

# Data required

The data used for the analysis of the maximum static density difference of temperature inverted profiles is from [ECCOv4r4](https://ecco-group.org/products-ECCO-V4r4.htm)[^1] project and the [GOSHIP easy ocean data product](https://cchdo.ucsd.edu/products/goship-easyocean)[^2].

## ECCOv4r4

From ECCOv4r4 the daily averaged potential temperature and practical salinity on the $0.5^{\circ}$ latitude-longitude grid is used.

## GOSHIP

From GOSHIP the reported profile data is used.

[^1]: G. Forget et al. “ECCO version 4: An integrated framework for non-linear inverse modeling and global ocean state estimation”. In: Geoscientific Model Development 8.10 (Oct. 2015), pp. 3071–3104. issn: 19919603. doi: 10.5194/gmd-8-3071-2015.

[^2]: Katsuro Katsumata et al. “GO-SHIP Easy Ocean: Gridded ship-based hydrographic section of temperature, salinity, and dissolved oxygen”. In: Scientific Data 9.1 (Dec. 2022). issn: 20524463. doi: 10.1038/s41597-022-01212-w.
