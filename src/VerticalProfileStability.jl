module VerticalProfileStability

using Rasters, OceanRasterConversions, GibbsSeaWater, Glob, Reexport

@reexport using Rasters, OceanRasterConversions, GibbsSeaWater, Glob, Dates

export
    ECCO_datadir,
    sim_datadir,
    plotdir,
    srcdir

const ECCO_datadir = joinpath(@__DIR__, "../data/observations/ECCO_daily_mean_TS/")
const sim_datadir = joinpath(@__DIR__, "../data/simulations/")
const plotdir = joinpath(@__DIR__, "../plots")
const srcdir = joinpath(@__DIR__)

include("maxdensitydiff.jl")
include("watercolmodel.jl")

@reexport using .VerticalProfileStability.MaximumDensityDifference
@reexport using .VerticalProfileStability.WaterColumnModel

end
