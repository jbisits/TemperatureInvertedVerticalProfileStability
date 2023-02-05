using CairoMakie
using .VerticalProfileStability

# Trying to write `convert_arguments` for a `Raster`

function Makie.convert_arguments(rs::Raster)

    lon, lat = collect(lookup(rs, X)), collect(lookup(rs, Y))
    plot_var = Matrix(rs[:, :, 1])

    return (lon, lat, plot_var)

end

output_series[Ti(1)][:Θᵤ]
convert_arguments(output_series[Ti(1)][:Θᵤ])

heatmap(output_series[Ti(1)][:Θᵤ])

## Expected output
heatmap(Matrix(output_series[Ti(1)][:Θᵤ][:, :, 1]))
