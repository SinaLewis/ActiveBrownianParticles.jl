using ActiveBrownianParticles, Random
using Plots
using Test

rng = Random.MersenneTwister()

posDist=rand(rng,Float64,(1000,2))
sum(posDist[:,1].^2)


plot(posDist[:,1],posDist[:,2],seriestype = :scatter)

vel = ActiveBrownianParticles()

@testset "ActiveBrownianParticles.jl" begin
    # Write your tests here.
end
