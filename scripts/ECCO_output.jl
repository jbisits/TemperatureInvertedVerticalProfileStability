## Saving files of variables from the monthly data
const ECCO_data_analysis = joinpath(@__DIR__, "../data/analysis/")

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
