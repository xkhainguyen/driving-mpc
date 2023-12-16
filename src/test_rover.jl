using Revise
using RobotZoo
import RobotDynamics as RD

model = RobotZoo.Rover(steer_rate_lim=0.5, accel_lim=5)
dmodel = RD.DiscretizedDynamics{RD.RK4}(model)
# model = RobotZoo.Rover()
n,m = size(model)

# Generate random state and control vector
# [x, y, θ, δ, v, a]
x = [0; 0; 0; 0; 20]  # initial condition
u = [2; 0.2]
dt = 0.005  # time step (s)
z = KnotPoint(x, u, 0, dt)
N = 100
xs = zeros(5,N)
xs[:,1] .= x
t = 0
# Evaluate the continuous dynamics and Jacobian
J = RD.DynamicsJacobian(model)
for i in 1:N-1
    # u = 3*rand(2) .- 0.5   
    xs[:,i+1]  = RD.discrete_dynamics(dmodel, xs[:,i] + 0.01*randn(5), u, t, dt)
end

display(xs) 

# Plotting
using Plots
# pyplot()
Plots.plot(xs[1,:], xs[2,:], xlabel="x position", ylabel="y position", markershape=:circle, reuse=false)