module ActiveBrownianParticles

#using Debugger
using Random, Plots

gr()

include("integration.jl")
include("equationOfMotion.jl")
include("updatePositions.jl")

# Define simulation params
N = 1000    # number of particles
boxSize = 25    # size of box

# Define random number generator
rng = MersenneTwister()

# Define starting points
# Randomly place in box
posDist = rand(rng, Float64, (N, 2)) .* boxSize
plot(posDist[:, 1], posDist[:, 2], seriestype=:scatter) # need display to show plot while debugging

updatePositions(10, 0.1, posDist, EOM, RK4step)

end
