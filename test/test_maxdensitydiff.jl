lon, lat, z = -180:5:180, -90:5:90, -400:100:0
timestamps = 1:2
Sₐ = Array{Float64}(undef, length(lon), length(lat), length(z))
fill!(Sₐ, 35)
Θ = Array{Float64}(undef, length(lon), length(lat), length(z))
Θ_vals = [-1.85, -1.0, -0.5, 0.25, 0.5]
for i ∈ eachindex(z)
    Θ[:, :, i] .= Θ_vals[i]
end

stack_vec = [RasterStack((Sₐ, Θ), (X(lon), Y(lat), Z(z)); name = (:SALT, :THETA))
             for i ∈ timestamps]
test_series = RasterSeries(stack_vec, Ti(timestamps))

ΔΘ_thres = 2.0
Δρ_max_series = series_max_Δρ(test_series, ΔΘ_thres; zdepth = -100.0)
