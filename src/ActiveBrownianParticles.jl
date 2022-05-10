using Random

module ActiveBrownianParticles
# Define simulation params
N = 1000    # number of particles
boxSize = 25    # size of box
T = 3.0     # temperature

# Define random number generator
rng = MersenneTwister()

# Define starting points
# Randomly place in box
posDist = rand(rng, Float64, (N, 2)) .* boxSize

# Define starting velocity
# Draw from normal distribution
velDist = randn(rng, Float64, (N, 2))

# determine scaling factor for temperature
KE = sum(velDist[:, 1] .^ 2 + velDist[:, 2] .^ 2)
λ = sqrt(N * T / KE)
# normalize velocity distriution to the temperature
velDist = λ .* velDist

######################



end
