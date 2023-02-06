using CairoMakie
using .VerticalProfileStability

# Trying to write `convert_arguments` for a `Raster`, this works!!
# Can write other recipes (scatter would be good if needed).

function Makie.convert_arguments(P::SurfaceLike, rs::Raster)

    lon, lat = collect(lookup(rs, X)), collect(lookup(rs, Y))
    plot_var = Matrix(rs[:, :, 1])

    return convert_arguments(P, lon, lat, plot_var)

end

heatmap(output_series[Ti(1)][:Θᵤ])
contourf(output_series[Ti(1)][:Θᵤ])

## Expected output
heatmap(Matrix(output_series[Ti(1)][:Θᵤ][:, :, 1]))
