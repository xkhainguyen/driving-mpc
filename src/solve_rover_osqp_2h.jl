"""OSQP with input constraint and partially state constraint
Once the horizon hits the end of the trajectory, copy the goal.
"""

# using Revise
using TrajectoryOptimization
import RobotZoo.Rover
using StaticArrays, LinearAlgebra
using Plots
using RobotDynamics
using Altro
const RD = RobotDynamics
using SparseArrays
using OSQP

include("utils.jl")

print("\nOSQP SIMULATION START!!!")

# Get reference trajectory
dt = 1/20   # 10 - 20Hz
N = 151
tf = (N-1)*dt

steer_rate_lim=0.3
accel_lim=3

model = Rover(ref=:rear)
n,m = size(model);
Nx,Nu = size(model);

x0 = @SVector zeros(n)
xf = @SVector [15, 15, π/4, 0, 0];  

# Objective function
Q = 1e-1*Diagonal(@SVector ones(n)) 
Qf = 1e3*Diagonal(@SVector ones(n))
R = 1e1*Diagonal(@SVector ones(m)) 
obj = LQRObjective(Q,R,Qf,xf,N);

# Create Empty ConstraintList
conSet = ConstraintList(n,m,N)

# Circle constraints
xc = SA[6, 10, 12]
yc = SA[5, 12, 5]
rc = SA[2, 2, 2]    #try 1
cir = CircleConstraint(n, xc, yc, rc)
add_constraint!(conSet, cir, 2:N-1)

# Control-state Bounds
u_bnd = [accel_lim; steer_rate_lim]
x_bnd = [Inf;Inf;Inf;40;π/3]
bnd = BoundConstraint(n,m, x_min=-x_bnd, x_max=x_bnd, u_min=-u_bnd, u_max=u_bnd)
add_constraint!(conSet, bnd, 1:N-1)
# Goal Constraint
goal = GoalConstraint(xf,1:2)
add_constraint!(conSet, goal, N)

prob = Problem(model, obj, x0, tf, xf=xf, constraints=conSet)

u0 = @SVector fill(0.0,m)
U0 = [u0 for k = 1:N-1]
initial_controls!(prob, U0)
rollout!(prob);

solver = ALTROSolver(prob)
set_options!(solver, show_summary=false)
initial_controls!(solver, U0)
solve!(solver);

Usln = hcat(Vector.(controls(solver))...)
Xsln = hcat(Vector.(states(solver))...)

# Rover simulator
xs = zeros(5,N)
x0_off = [0;1;0;0;0]
xs[:,1] .= x0 + x0_off
dd = zeros(5,N)

# dmodel2 = RD.DiscretizedDynamics{RD.RK4}(model2)
for k in 1:N-1
    dd[:,k] = 0*randn(5)
    # xs[:,k+1] .= simulator(model, xs[:,k]+dd[:,k], Usln[:,k], u_bnd, dt)
    xs[:,k+1] .= DT_non_model(model, xs[:,k]+dd[:,k], Usln[:,k], dt)
end

Q = sparse(1.0*I(Nx))
R = sparse(1.0*I(Nu))
Qf = sparse(1.0*I(Nx))
  
# Get cost-to-go from Riccati
K, P = riccati_recursion(model, Xn, Un)

thist = Array(range(0,dt*(N-1), step=dt))

# Initialize OSQP
Nh = 30
Nh′ = Nh

Un = copy(Usln)
Xn = copy(Xsln)  # ref state

# copy the end point to satisfy horizon
for k = 1:Nh
    Un = hcat(Un, Un[:,end])
    Xn = hcat(Xn, Xn[:,end])
    P = cat(P, P[:,:,end], dims=(3,3))
end

xhist = zeros(n,N)
xhist[:,1] .= xs[:,1]
uhist = zeros(m,N-1)

for k = 1:(N-1)
    println(k)

    ubar0 = mpc_controller(k, xhist[:,k], Nh, model, P, Xn, Un, dt, x_bnd, u_bnd)

    uhist[:,k] .= Un[:,k] + ubar0
    
    # xhist[:,k+1] .= DT_non_model(model, xhist[:,k], uhist[:,k], dt)

    Ak1, Bk1 = DT_lin_model(model,Xn[:,k],Un[:,k],dt)
    xhist[:,k+1] .= Xn[:,k+1] + Ak1*(xhist[:,k] - Xn[:,k]) + Bk1*ubar0
end

# Plot solution
plot_scene("OSQP3.5")

print("OSQP SIMULATION FINISHED!!!")

##
include("visualization.jl")
vis = rover_vis_init()
##
rover_vis_run(vis, Xsln, xhist, Nh, xc, yc, rc, xf, visextralines=true)
