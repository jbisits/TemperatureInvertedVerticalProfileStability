module VerticalProfileStability

using Rasters, OceanRasterConversions, GibbsSeaWater, Glob, CairoMakie, Reexport

@reexport using Rasters, OceanRasterConversions, GibbsSeaWater, Glob, Dates, CairoMakie

export
    ECCO_DATADIR,
    SIM_DATADIR,
    PLOTDIR,
    SRCDIR

const ECCO_DATADIR = joinpath(@__DIR__, "../data/observations/ECCO_daily_mean_TS/")
const SIM_DATADIR = joinpath(@__DIR__, "../data/simulations/")
const PLOTDIR = joinpath(@__DIR__, "../plots")
const SRCDIR = joinpath(@__DIR__)

include("maxdensitydiff.jl")
include("watercolmodel.jl")
include("makierasterplot.jl")

@reexport using .VerticalProfileStability.MaximumDensityDifference
@reexport using .VerticalProfileStability.WaterColumnModel

end
