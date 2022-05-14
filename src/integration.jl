"""
Perform Runge-Kutta with the given function and initial condition
Takes the time, the array of particle positions, the step size, and the function of the
equation of motion
"""
function RK4step(t0, y0, stepSize::Float64, f::Function)
    # given an initial value problem
    # dy/dt = f(t,y), y(t_0) = y_0

    # pick a step size
    h = stepSize

    # define the slope at some points
    k1 = zeros(length(y0[:, 1]), 2)
    k2 = zeros(length(y0[:, 1]), 2)
    k3 = zeros(length(y0[:, 1]), 2)
    k4 = zeros(length(y0[:, 1]), 2)

    for i in 1:length(y0[:, 1])
        k1[i, :] = f(t0, y0, i) # initial point
        k2[i, :] = f(t0 + h / 2, y0 .+ h * k1[i] / 2, i) # middle point
        k3[i, :] = f(t0 + h / 2, y0 .+ h * k2[i] / 2, i) # middle point again
        k4[i, :] = f(t0 + h, y0 .+ h * k3[i], i)     # end point
    end

    # update positions
    y1 = y0 .+ (h / 6) .* (k1 .+ (2 .* k2) .+ (2 .* k3) .+ k4)
    t1 = t0 + h
    return t1, y1
end