## Saving files of variables from the monthly data
const ECCO_data_analysis = joinpath(@__DIR__, "../data/analysis/")
using Glob, Dates
using .VerticalProfileStability

ECCO_data = glob("*.nc", ECCO_datadir)
timestamps = Date(2007, 01, 01):Day(1):Date(2007, 12, 31)
series = RasterSeries(ECCO_data, Ti(timestamps); child = RasterStack)

## ΔΘ_thres = 0.5
ΔΘ_thres = 0.5
mkdir(joinpath(ECCO_data_analysis, "output_$ΔΘ_thres"))
savepath = joinpath(ECCO_data_analysis, "output_$ΔΘ_thres")
series_max_Δρ(series, ΔΘ_thres, savepath)

## ΔΘ_thres = 1.0
ΔΘ_thres = 1.0
mkdir(joinpath(ECCO_data_analysis, "output_$ΔΘ_thres"))
savepath = joinpath(ECCO_data_analysis, "output_$ΔΘ_thres")
series_max_Δρ(series, ΔΘ_thres, savepath)

## ΔΘ_thres = 2.0
ΔΘ_thres = 2.0
mkdir(joinpath(ECCO_data_analysis, "output_$ΔΘ_thres"))
savepath = joinpath(ECCO_data_analysis, "output_$ΔΘ_thres")
series_max_Δρ(series, ΔΘ_thres, savepath)

## 1.0 ≤ ΔΘ_thres < 2.5
ΔΘ_thres = [1.0, 2.5]
mkdir(joinpath(ECCO_data_analysis, "output_$(ΔΘ_thres[1])lower_$(ΔΘ_thres[2])upper"))
savepath = joinpath(ECCO_data_analysis, "output_$(ΔΘ_thres[1])lower_$(ΔΘ_thres[2])upper")
series_max_Δρ(series, ΔΘ_thres, savepath)
