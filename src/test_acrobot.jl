using Revise
using RobotZoo
using RobotDynamics
using ForwardDiff
using FiniteDiff

model = RobotZoo.Cartpole()
n,m = size(model)

# Generate random state and control vector
x,u = rand(model)
dt = 0.1  # time step (s)
z = KnotPoint(x,u,0,dt)

# Evaluate ythe continuous dynamics and Jacobian
ẋ = RobotDynamics.dynamics(model, x, u)
∇f = RobotDynamics.DynamicsJacobian(model)
RobotDynamics.jacobian!(model, ∇f, ẋ, z)

# Evaluate the discrete dynamics and Jacobian
x′ = discrete_dynamics(RK3, model, z)
discrete_jacobian!(RK3, ∇f, model, z)