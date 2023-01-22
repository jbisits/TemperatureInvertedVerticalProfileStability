module VerticalProfileStability

using Rasters, OceanRasterConversions, GibbsSeaWater, Reexport

@reexport using Rasters, OceanRasterConversions, GibbsSeaWater

include("max_density_diff.jl")

end
