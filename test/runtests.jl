using Test

examples_dir = joinpath(@__DIR__, "..", "examples")
@assert ispath(examples_dir)
example_filenames = readdir(examples_dir)
@assert "ft01_poisson.jl" in example_filenames

@testset "$(filename)" for filename in example_filenames
    path = joinpath(examples_dir, filename)
    @test (include(path) ;true)
end
