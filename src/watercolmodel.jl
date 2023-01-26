module WaterColumnModel

export
    run_model,
    run_model_two_layer,
    order_files,
    order_mat,
    vec2mat,
    save_timeseries!,
    save_Δρ_timeseries!

using Oceananigans, SeawaterPolynomials.TEOS10, JLD2, Glob, GibbsSeaWater
using Oceananigans.Units: seconds, minutes, hours, days
using Oceananigans: ∂z_b

############################################################################################
## Water column model
############################################################################################
"""
    function run_model(; params)
Run the model which is an ensemble of water column experiments, each with a
different initial condition for salinity in the mixed layer.
The salinity initial conditions set in the mixed layer which are determined by
        `S_ics = range(Sᵤ, Sᵣ, length = num_ics)`
here `Sᵣ` is the upper limit for salinity initial condition in the mixed layer.
The default values are what I will use for my current simulations
but setting things up like this means they can be changed easily.
For information on how the density differences are calculated see `model_setup.md`.
"""
function run_model(;
                    Tᵤ = -1.85,
                    Tₗ = 0.5,
                    Sᵤ = 34.5,
                    Sₗ = 34.7,
                    Sₘ = 34.75,
                    Sᵣ = 34.57,
               num_ics = 10,
                    Nz = 100, # number of levels
                    Lz = 500, # ovearll depth
           ref_density = 1028.18,
          convective_κ = 1.0,
          background_κ = 1e-5,
            convective = true,
                    Δt = 5,   # minutes
              savepath = "Output",
            sim_length = 60,  # in days
             save_freq = 30   # in minutes
        )

    # Save the parameters
    saved_params = joinpath(savepath, "model_params.jld2")
    jldsave(saved_params;  Tᵤ, Tₗ, Sᵤ, Sₗ, Sₘ, Sᵣ, num_ics, sim_length, save_freq)

    # Grid
    grid = RectilinearGrid(size = Nz, z = (-Lz, 0), topology=(Flat, Flat, Bounded))

    # Buoyancy, using TEOS10
    EOS = TEOS10EquationOfState(; reference_density = ref_density)
    buoyancy = SeawaterBuoyancy(equation_of_state = EOS)

    ## Turbulence closure
    closure = convective==true ? ConvectiveAdjustmentVerticalDiffusivity(
                                 convective_κz = convective_κ,
                                 background_κz = background_κ) :
                                 ScalarDiffusivity(κ = background_κ)

    # Set temperature initial condition
    T₀ = Array{Float64}(undef, size(grid))
    Tᵤ_array = fill(Tᵤ, 20) # want to unhard code these bits eventually
    # This adds a temperature gradient to avoid spurios convective mixing in the mixed layer
    Tₗ_array = fill(Tₗ, 80)
    #Tᵤ_array = reverse(range(-1.87, Tᵤ, length = 20))
    T₀[:, :, :] = vcat(Tₗ_array, Tᵤ_array)

    # Set the salinity initial condition
    S₀ = Array{Float64}(undef, size(grid))
    Sᵤ_array = fill(Sᵤ, 20)
    Sₗ_array = range(Sₘ, Sₗ, length = 80)
    S₀[:, :, :] = vcat(Sₗ_array, Sᵤ_array)
    # Set the salinity initial conditions for increasing salinity in the mixed layer
    Sᵢ = reshape(S₀, :)
    S_ics = repeat(Sᵢ, inner = (1, 1), outer = (1, num_ics))
    S_inc = range(Sᵤ, Sᵣ, length = num_ics)
    S_ics[81:end, :] .= repeat(S_inc, inner = (1, 20), outer = (1, 1))'

    savefile = joinpath(savepath, "T_and_S")
    for i ∈ 1:num_ics

        @info "Simulation number $i"
        model = NonhydrostaticModel(grid = grid,
                                    tracers = (:T, :S),
                                    coriolis = FPlane(f = 1e-4),
                                    buoyancy = buoyancy,
                                    closure = closure)

        S₀ = Array{Float64}(undef, size(model.grid))
        S₀[:, :, :] = S_ics[:, i]

        set!(model, T = T₀, S = S₀)

        simulation = Simulation(model, Δt = Δt * minutes, stop_time = sim_length * days)

        outputs = convective==true ? (T = model.tracers.T, S = model.tracers.S,
                                      κ = save_diffusivity, ∂z_b = save_∂z_b) :
                                     (T = model.tracers.T, S = model.tracers.S)
        simulation.output_writers[:outputs] = JLD2OutputWriter(model, outputs,
                                               filename = savefile*string(i)*".jld2",
                                               schedule = TimeInterval(save_freq * minutes))

        run!(simulation)
    end
end
"""
A two layer version of the same model with changes where appropriate.
"""
function run_model_two_layer(;
                            Tᵤ = -1.85,
                            Tₗ = 0.5,
                            Sᵤ = 34.5,
                            Sₗ = 34.7,
                            Sᵣ = 34.578,
                       num_ics = 11,
                            Nz = 2,         # number of levels
                            Lz = 200,       # ovearll depth
                   ref_density = 1028.18,
                  convective_κ = 1.0,
                  background_κ = 1e-5,
                    convective = true,
                            Δt = 5,         # minutes
                      savepath = "Output",
                    sim_length = 60,        # in days
                     save_freq = 30         # in minutes
                        )

    # Save the parameters
    saved_params = joinpath(savepath, "model_params.jld2")
    jldsave(saved_params;  Tᵤ, Tₗ, Sᵤ, Sₗ, Sᵣ, num_ics, sim_length, save_freq)

    # Grid
    grid = RectilinearGrid(size = Nz, z = (-Lz, 0), topology=(Flat, Flat, Bounded))

    # Buoyancy, using TEOS10
    EOS = TEOS10EquationOfState(; reference_density = ref_density)
    buoyancy = SeawaterBuoyancy(equation_of_state = EOS)

    ## Turbulence closure
    closure = convective==true ? ConvectiveAdjustmentVerticalDiffusivity(
                                                    convective_κz = convective_κ,
                                                    background_κz = background_κ) :
                                                    ScalarDiffusivity(κ = background_κ)

    # Set temperature initial condition
    T₀ = Array{Float64}(undef, size(grid))
    Tᵤ_array = fill(Tᵤ, 1)
    Tₗ_array = fill(Tₗ, 1)
    T₀[:, :, :] = vcat(Tₗ_array, Tᵤ_array)

    # Set the salinity initial condition
    S₀ = Array{Float64}(undef, size(grid))
    Sᵤ_array = fill(Sᵤ, 1)
    Sₗ_array = fill(Sₗ, 1)
    S₀[:, :, :] = vcat(Sₗ_array, Sᵤ_array)
    # Set the salinity initial conditions for increasing salinity in the mixed layer
    Sᵢ = reshape(S₀, :)
    S_ics = repeat(Sᵢ, inner = (1, 1), outer = (1, num_ics))
    S_inc = range(Sᵤ, Sᵣ, length = num_ics)
    S_ics[2, :] .= repeat(S_inc, inner = (1, 1), outer = (1, 1))

    # Save the salinity mixed layer initial conditions
    jldopen(saved_params, "a+") do file
        file["T_ics"] = reshape(T₀, :)
        file["S_ics"] = S_ics
    end

    savefile = joinpath(savepath, "T_and_S")
    for i ∈ 1:num_ics

    @info "Simulation number $i"
    model = NonhydrostaticModel(grid = grid,
                        tracers = (:T, :S),
                        coriolis = FPlane(f = 1e-4),
                        buoyancy = buoyancy,
                        closure = closure)

    S₀ = Array{Float64}(undef, size(model.grid))
    S₀[:, :, :] = S_ics[:, i]

    set!(model, T = T₀, S = S₀)

    simulation = Simulation(model, Δt = Δt * minutes, stop_time = sim_length * days)

    outputs = convective==true ? (T = model.tracers.T, S = model.tracers.S,
                                  κ = save_diffusivity, ∂z_b = save_∂z_b) :
                                 (T = model.tracers.T, S = model.tracers.S)
    simulation.output_writers[:outputs] = JLD2OutputWriter(model, outputs,
                                filename = savefile*string(i)*".jld2",
                                schedule = TimeInterval(save_freq * minutes))

    run!(simulation)
    end
end
"""
    save_∂z_b(model)
Save the buoyancy gradient that is calculated during the simulations.
The `ConvectiveAdjustmentVerticalDiffusivity` closure computes the buoyancy gradient using
`∂z_b` and applies the convective diffusivity if `∂z_b < 0` and the background diffusivity
if `∂z_b ≥ 0`.
"""
function save_∂z_b(model)

    buoyancy_gradient = Array{Float64}(undef, model.grid.Nz)

    for i ∈ 1:model.grid.Nz
        buoyancy_gradient[i] = ∂z_b(1, 1, i, model.grid, model.buoyancy, model.tracers)
    end

    return buoyancy_gradient
end
"""
    save_diffusivity(model)
Save the diffusivity that is applied to a region using the method in the
`ConvectiveAdjustmentVerticalDiffusivity` closure. The code obtained to do this is from
`Oceananigans.jl` though I do not use the `@kernel` macro as I am not on a `GPU`.
"""
function save_diffusivity(model)

    is_stableᶜᶜᶠ(i, j, k, grid, tracers, buoyancy)=∂z_b(i, j, k, grid, buoyancy, tracers)>=0
    diffusivities = Array{Float64}(undef, model.grid.Nz)

    for i ∈ 1:model.grid.Nz
        stable_cell = is_stableᶜᶜᶠ(1, 1, i, model.grid, model.tracers, model.buoyancy)

        diffusivities[i] = ifelse(stable_cell,
                                  model.closure.background_κz,
                                  model.closure.convective_κz)
    end

    return diffusivities
end

############################################################################################
## Functions for organising and analysing model output
############################################################################################

"""
    function Δρ(data_files::Vector{String}, ΔΘ_thres::Float64, p_ref::Int64)
Find the greatest denstiy difference over the water column from the ouput of a
one dimensional model. This method uses a reference pressure, `p_ref`, instead of
midpoint or lower pressure for the density difference calculations. The simulations I will
run as part of this will use reference pressure so only have this method here.
"""
function Δρ(data_files::Vector{String}, ΔΘ_thres::Float64, p_ref::Int64)

    T_test = FieldTimeSeries(data_files[1], "T")
    t = T_test.times ./ days
    z = znodes(T_test)
    all_levels = reverse(z)
    num_steps = length(t)

    Δρ_s = Array{Float64}(undef, num_steps, length(data_files))
    Δρ_c = Array{Float64}(undef, num_steps, length(data_files))

    for (i, file) ∈ enumerate(data_files)

        println("$file")
        T, S = FieldTimeSeries(file, "T"), FieldTimeSeries(file, "S")

        #T and S time series for each initial condition
        T_steps = [interior(T[time], 1, 1, :) for time ∈ 1:length(t)]
        S_steps = [interior(S[time], 1, 1, :) for time ∈ 1:length(t)]

        for step ∈ 1:num_steps

            # The arrays for depth, pressure have been reversed as it is
            # conceptually easier to work with so T and S have to be as well
            T_profile = reverse(T_steps[step])
            S_profile = reverse(S_steps[step])
            temp_Δρ_s = []
            temp_Δρ_c = []
            for upper_level ∈ eachindex(all_levels[1:end-2])

                Θᵤ = T_profile[upper_level]
                Sᵤ = S_profile[upper_level]

                for lower_level ∈ (upper_level + 1):length(all_levels)

                    Θₗ = T_profile[lower_level]
                    ΔΘ = Θᵤ - Θₗ
                    if  abs(ΔΘ) > ΔΘ_thres
                        Sₗ = S_profile[lower_level]
                        ΔS = Sᵤ - Sₗ
                        ρᵤ = gsw_rho(Sᵤ, Θᵤ, p_ref)
                        ρₗ = gsw_rho(Sₗ, Θₗ, p_ref)
                        Δρ_s_level = ρᵤ - ρₗ
                        push!(temp_Δρ_s, Δρ_s_level)
                        α = gsw_alpha(Sₗ, Θₗ, p_ref)
                        β = gsw_beta(Sₗ, Θₗ, p_ref)
                        Δρ_c_level = ρₗ * (β * ΔS - α * ΔΘ)
                        push!(temp_Δρ_c, Δρ_c_level)
                    end

                end

            end

            if !isempty(temp_Δρ_s)
                Δρ_s[step, i] = maximum(temp_Δρ_s)
                Δρ_c[step, i] = maximum(temp_Δρ_c)
            end

        end

    end

    return Dict("Δρ_s" => Δρ_s,
                "Δρ_c" => Δρ_c)
end

"""
    function order_files(files::Vector{String})
Re order the files so that they are increasing numerically not alphabetically.
"""
function order_files(files::Vector{String})

    num_files = length(files)
    diff_from_10 = num_files - 10
    first_9 = vcat(1, (3+diff_from_10):num_files)
    ordered_files = vcat(files[first_9], files[2:(2+diff_from_10)])

    return ordered_files
end

"""
    function order_mat(data::Matrix)
Reorder the columns of the data generated by the Δρ function so that
it is in order of most stable to most unstable.
"""
function order_mat(data::Matrix)

    first_9 = vcat(1, 3:10)
    ordered_matrix = hcat(data[:, first_9], data[:, 2])

    return ordered_matrix

end

"""
    function vec2mat(data::Vector)
Turn the output data which is saved as an array of vectors into a matrix.
"""
function vec2mat(data)

    mat = [data[1] data[2]]
    length = 3:lastindex(data)
    for i ∈ length
        mat = hcat(mat, data[i])
    end

    return mat

end

"""
    function extract_t(data_files::Vector{String})
Extract the time (`t`) for the time series from the output files.
"""
function extract_t(data_files::Vector{String}; time_length = days)
    file1 = jldopen(data_files[1])
    iterations = parse.(Int, keys(file1["timeseries/t"]))
    t = [file1["timeseries/t/$i"] for i ∈ iterations] ./ time_length
    close(file1)
    return t
end

"""
    function extract_T_ts(data_files::Vector{String})
Extract the temperature (`T`) time series from the output files.
"""
function extract_T_ts(data_files::Vector{String})

    file1 = jldopen(data_files[1])
    iterations = parse.(Int, keys(file1["timeseries/t"]))
    num_levels = file1["grid/Nz"]
    close(file1)

    T_ts = Array{Float64}(undef, num_levels, length(iterations), length(data_files))
    ic = 1
    for file ∈ data_files
        T = FieldTimeSeries(file, "T")
        T_ts[:, :, ic] = vec2mat([interior(T[i], 1, 1, :) for i ∈ eachindex(iterations)])
        ic += 1
    end

    return T_ts
end

"""
    function extract_S_ts(data_files::Vector{String})
Extract the salinity (`S`) time series from the output files.
"""
function extract_S_ts(data_files::Vector{String})

    file1 = jldopen(data_files[1])
    iterations = parse.(Int, keys(file1["timeseries/t"]))
    num_levels = file1["grid/Nz"]
    close(file1)

    S_ts = Array{Float64}(undef, num_levels, length(iterations), length(data_files))
    ic = 1
    for file ∈ data_files
        S = FieldTimeSeries(file, "S")
        S_ts[:, :, ic] = vec2mat([interior(S[i], 1, 1, :) for i ∈ eachindex(iterations)])
        ic += 1
    end

    return S_ts
end

"""
    function extract_κ_ts(data_files::Vector{String})
Extract the diffusivity (`κ`) time series from the output files.
"""
function extract_κ_ts(data_files::Vector{String})

    file1 = jldopen(data_files[1])
    iterations = parse.(Int, keys(file1["timeseries/t"]))
    num_levels = file1["grid/Nz"]
    close(file1)

    κ_ts = Array{Float64}(undef, num_levels, length(iterations), length(data_files))
    ic = 1
    for file ∈ data_files
        open_file = jldopen(file)
        κ = [open_file["timeseries/κ/$i"] for i ∈ iterations]
        close(open_file)
        κ_ts[:, :, ic] = vec2mat(κ)
        ic += 1
    end

    return κ_ts
end

"""
    function extract_∂z_b_ts(data_files::Vector{String})
Extract the buoyancy gradient (`∂z_b`) time series from the output files.
"""
function extract_∂z_b_ts(data_files::Vector{String})

    file1 = jldopen(data_files[1])
    iterations = parse.(Int, keys(file1["timeseries/t"]))
    num_levels = file1["grid/Nz"]
    close(file1)

    ∂z_b_ts = Array{Float64}(undef, num_levels, length(iterations), length(data_files))
    ic = 1
    for file ∈ data_files
        open_file = jldopen(file)
        ∂z_b_temp = [open_file["timeseries/∂z_b/$i"] for i ∈ iterations]
        close(open_file)
        ∂z_b_ts[:, :, ic] = vec2mat(∂z_b_temp)
        ic += 1
    end

    return ∂z_b_ts
end

"""
    function save_timeseries!(filename::String, data_files::Vector{String};
                                ΔΘ_thres = 0.5, p_ref = 0)
Save the time series of `T`, `S`, `κ`, `∂z_b` and Δρ (the maximum density difference where the
temperature difference between the upper and lower level (ΔΘ) is greater than `ΔΘ_thres`)
from the output data. All `GibbsSeaWater` computations use reference pressure p_ref.
"""
function save_timeseries!(filename::String, data_files::Vector{String};
                            ΔΘ_thres = 0.5, p_ref = 0)

    @info "Time"
    t = extract_t(data_files)
    @info "Temperature"
    T_ts = extract_T_ts(data_files)
    @info "Salinity"
    S_ts = extract_S_ts(data_files)
    @info "Diffusivity"
    κ_ts = extract_κ_ts(data_files)
    @info "Buoyancy gradient"
    ∂z_b_ts = extract_∂z_b_ts(data_files)
    @info "Maximum density difference"
    Δρ_ts = Δρ(data_files, ΔΘ_thres, p_ref)

    jldsave(filename; t, T_ts, S_ts, κ_ts, ∂z_b_ts, Δρ_ts)
    return nothing
end
"""
    function save_Δρ_timeseries!(filename::String, data_files::Vector{String};
                                 ΔΘ_thres = 0.5, p_ref = 0)
Only save the density difference timeseries. This is to more quickly look at different
values for for `ΔΘ_thres`.
"""
function save_Δρ_timeseries!(filename::String, data_files::Vector{String};
                             ΔΘ_thres = 0.5, p_ref = 0)

    @info "Maximum density difference"
    Δρ_ts = Δρ(data_files, ΔΘ_thres, p_ref)

    jldsave(filename; Δρ_ts)
    return nothing
end

end # module
