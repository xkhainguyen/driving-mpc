using Revise
using RobotZoo
using RobotDynamics

model = RobotZoo.Rover(steer_lim=0.5, accel_lim=5)
# model = RobotZoo.Rover()
n,m = size(model)

# Generate random state and control vector
# [x, y, θ, δ, v, a]
x = [0; 0; 0; 0; 20; 0]  # initial condition
u = [2; 0.2]
dt = 0.005  # time step (s)
z = KnotPoint(x, u, 0, dt)
N = 20
xs = zeros(6,N)
xs[:,1] .= x
# Evaluate the continuous dynamics and Jacobian
J = RobotDynamics.DynamicsJacobian(model)
for i in 1:N-1
    # u = 3*rand(2) .- 0.5   
    ẋ = RobotDynamics.dynamics(model, xs[:,i], u)
    # RobotDynamics.jacobian!(model, J, ẋ, z)
    xs[:,i+1] = ẋ*dt + xs[:,i];
end


# model_d = RobotDynamics.DiscretizedDynamics(model, RobotDynamics.RK3)

# Evaluate the discrete dynamics and Jacobian
# x′ = RobotDynamics.discrete_dynamics(model, z) 
# discrete_jacobian!(RK3, ∇f, model, z)

display(xs) 

# Plotting
using Plots
Plots.plot(xs[1,:], xs[2,:], xlabel="x position", ylabel="y position", markershape=:circle)