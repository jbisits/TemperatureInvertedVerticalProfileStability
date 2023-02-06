using .VerticalProfileStability
# Trying to write `convert_arguments` for a `Raster`, this works!!
# Can write other recipes (scatter would be good if needed).

## ::SurfaceLike
function Makie.convert_arguments(P::SurfaceLike, rs::Raster)

    lon, lat = collect(lookup(rs, X)), collect(lookup(rs, Y))
    plot_var = Matrix(rs[:, :, 1])

    return convert_arguments(P, lon, lat, plot_var)

end

## ::PointBased
function Makie.convert_arguments(P::PointBased, rs_1::Raster, rs_2::Raster)

    var_1_vec = collect(reshape(rs_1, :))
    var_2_vec = collect(reshape(rs_2, :))

    return convert_arguments(P, var_1_vec, var_2_vec)

end

## Check the `SurfaceLike` and `PointBased` plots,
# `output_series` is from `../data/analysis/ECCO_output_analysis.jl`

# heatmap(output_series[Ti(1)][:Θᵤ])
# contourf(output_series[Ti(1)][:Θᵤ])
# scatter(output_series[Ti(1)][:Θₗ], output_series[Ti(1)][:Δρ_cab])

## Example of plotting for a series
# fig = Figure()
# ax = Axis(fig[1, 1])
# for s ∈ output_series[Ti(1:2)]
#     scatter!(s[:Θₗ], s[:Δρ_cab])
# end
# fig
