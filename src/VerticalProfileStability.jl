module VerticalProfileStability

using Rasters, OceanRasterConversions, GibbsSeaWater, Glob, CairoMakie, JLD2, Reexport

@reexport using Rasters,
                OceanRasterConversions,
                GibbsSeaWater,
                Glob,
                Dates,
                CairoMakie,
                NCDatasets,
                JLD2

export
    ECCO_DATADIR,
    EN4_DATADIR,
    SIM_DATADIR,
    GOSHIP_DATADIR,
    PLOTDIR,
    SRCDIR

const ECCO_DATADIR = joinpath(@__DIR__, "../data/observations/ECCO_daily_mean_TS/")
const EN4_DATADIR = joinpath(@__DIR__, "../data/observations/EN4/")
const GOSHIP_DATADIR = joinpath(@__DIR__, "../data/observations/goship/mat_files")
const SIM_DATADIR = joinpath(@__DIR__, "../data/simulations/")
const PLOTDIR = joinpath(@__DIR__, "../plots")
const SRCDIR = joinpath(@__DIR__)

include("maxdensitydiff.jl")
include("watercolmodel.jl")
include("makierasterplot.jl")

@reexport using .VerticalProfileStability.MaximumDensityDifference
@reexport using .VerticalProfileStability.WaterColumnModel

end
