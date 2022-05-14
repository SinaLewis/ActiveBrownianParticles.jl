using Distributions

"""
Define the equation of motion
From the Redner, Hagan, Baskaran, Phys. Rev. Lett. 2013 paper
Takes the time, array of particle positions, and the particle whose EOM we're solving
"""
function EOM(t, y, i)
    σ = 3.0 * 10e-10  # units of meters
    D = 1.0 * 10e-9   # order of magnitude correct, in units of m^2/s
    Dr = 3 * D / σ^2  # valid in low Reynolds number limit
    β = 1.0   # in units of kB*T, is about 8.6*10E-5 eV

    # Repulsive interparticle interaction
    # Uses a WCA potential
    # V = 4ϵ[(σ/r)^12 - (σ/r)^6] + ϵ for r < 2^(1//6)
    # F = -dV/dr
    # Use dimensionless length l = r/σ, only applies for l < 2^(1/6)/σ
    # Use ϵ=kB*T
    Fex = (0, 0)
    for j in 1:length(y[:, 1])
        l = sqrt(sum((y[i, :] .- y[j, :]) .^ 2)) / σ
        if j != i && l < 2^(1 / 6)
            # magnitude and vector pointing from j to i
            r = ((y[i, 1] - y[j, 1], y[i, 2] - y[j, 2]) ./ l)
            Fex = Fex .+ (4 * D) * (12 * l^(-13) - 6 * l^(-7)) .* r
        end
    end

    # random force
    mu = 0    #The mean of the Normal
    sigma = sqrt(2 * D) #The standard deviation of the Normal
    lb = -3 * sigma    #The truncation lower bound
    ub = 3 * sigma    #The truncation upper bound
    d = Truncated(Normal(mu, sigma), lb, ub)
    η = sqrt(2 * D) .* rand(d, 2)

    # self-propulsion force
    mu = 0    #The mean of the Normal
    sigma = sqrt(2 * Dr) #The standard deviation of the Normal
    lb = -3 * sigma    #The truncation lower bound
    ub = 3 * sigma    #The truncation upper bound
    d = Truncated(Normal(mu, sigma), lb, ub)
    θ = sqrt(2 * Dr) .* rand(d)

    Fp = (D * β) .* (cos(θ), sin(θ))

    ###### TOTAL 'FORCE' ##########
    F = Fex .+ Fp .+ η
end