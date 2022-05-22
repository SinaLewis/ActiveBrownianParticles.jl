"""
Update the positions of particles using the given equation of motion.
Updates the positions y0, n times, using the equation of motion f and timestep δt.
"""

function updatePositions(n::Int64, δt, y0, f::Function, solver::Function)

    for k in 1:n #number of time steps
        t0, y0 = solver(k, y0, δt, f)
        display(plot(y0[:, 1], y0[:, 2], seriestype=:scatter))
    end

    return y0
end
