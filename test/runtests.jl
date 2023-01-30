using Test, Rasters, OceanRasterConversions, Glob, JLD2
using .VerticalProfileStability

# @testset "Maximum density difference" begin

# end

include("test_wcmodel.jl")
@testset "Water column model" begin

    @test isequal(T_test[1, 1], T_test[1, 2])
    @test isequal(T_test[2, 1], T_test[2, 2])
    @test isequal(S_test[1, 2], S_test[1, 2])
    @test isequal(S_test[2, 2], S_test[2 ,2])

end
