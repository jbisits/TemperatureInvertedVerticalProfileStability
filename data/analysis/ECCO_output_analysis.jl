using CairoMakie
using .VerticalProfileStability

timestamps = Date(2007, 01, 01):Day(1):Date(2007, 12, 31)
ΔΘ_thres = [0.5, 1.0, 2.0]
output_path = joinpath(@__DIR__, "output_$(ΔΘ_thres[3])")
output_files = glob("*.nc", output_path)
output_series = RasterSeries(output_files, Ti(timestamps); child = RasterStack)
