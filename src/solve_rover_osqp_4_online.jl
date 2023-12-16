"""OSQP with input constraint and obstacle constraint
Once the horizon hits the end of the trajectory, copy the goal. Try warm-start
Linearize around previous solution
"""

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
include("my_mpc.jl")

printstyled("\nOSQP SIMULATION START!!!\n", bold=true, color=:green)

# Get reference trajectory
dt = 1/5   # 10 - 20Hz
N = 141
tf = (N-1)*dt

steer_rate_lim = 0.3
accel_lim = 3

model = Rover(ref=:rear)
n,m = size(model);
Nx, Nu = size(model);

x0 = @SVector zeros(Nx)
xf = @SVector [20, 20, π/4, 0, 0]

# Objective function
Q = 1e-1*Diagonal(@SVector ones(n)) 
Qf = 1e3*Diagonal(@SVector ones(n))
R = 1e1*Diagonal(@SVector ones(m)) 
obj = LQRObjective(Q,R,Qf,xf,N);

# Create Empty ConstraintList
conSet = ConstraintList(n,m,N)

# Circle constraints
xc = SA[14, 10, 16]
yc = SA[6, 10, 15]
rc = SA[1.5, 1.5, 1.5] 
rc_plan = SA[2.1, 2.1, 2.1]    #try 1

# Control-state Bounds
u_bnd = [accel_lim; steer_rate_lim]
x_bnd = [Inf;Inf;Inf;1.5;π/2.5]
bnd = BoundConstraint(n,m, x_min=-x_bnd, x_max=x_bnd, u_min=-u_bnd, u_max=u_bnd)
add_constraint!(conSet, bnd, 1:N-1)

prob = Problem(model, obj, x0, tf, xf=xf, constraints=conSet)

u0 = @SVector fill(0.0,m)
U0 = [u0 for k = 1:N-1]
initial_controls!(prob, U0)
rollout!(prob);

solver = ALTROSolver(prob)
set_options!(solver, show_summary=false)
initial_controls!(solver, U0)
solve!(solver);

Usln = Vector.(controls(solver))
Xsln = Vector.(states(solver))

Q = sparse(1.0*I(Nx))
R = sparse(1.0*I(Nu))
Qf = sparse(1.0*I(Nx))
  
Un = copy(Usln)
Xn = copy(Xsln)  # ref state

# Initialize OSQP
Nh = 40
# 10, 80

# copy the end point to satisfy horizon
for k = 1:Nh
    global Un = push!(Un, Un[end])
    global Xn = push!(Xn, Xn[end])
end

x_ws = zeros(Nh*(Nx+Nu))        # for warm start
obs_flag = [false for j=1:N]    # indicates when to see obstacles

xhist = zeros(n,N)
x0_off = [0;0.5;0;0;0]           # initial offset
xhist[:,1] = x0 + x0_off
uhist = zeros(m,N-1)

Xhorizon = [[zeros(Nx) for i=1:Nh] for j=1:N-1]
Uhorizon = [[zeros(Nu) for i=1:Nh] for j=1:N-1]
# store previous solution to feed online linear mpc, 
Xprev = [Xn[i] for i in 1:Nh+1]
Uprev = [Un[i] for i in 1:Nh]

u_bnd1 = 1.2*u_bnd

for k = 1:N-1
    println(k)

    Xhorizon[k], Uhorizon[k] = mpc_controller_online(k, xhist[:,k], xf, Nh, 
    model, Xprev, Uprev, dt, x_bnd, u_bnd1, warm_start=true)
    uhist[:,k] = Uhorizon[k][1] 

    xhist[:,k+1] = DT_non_model(model, xhist[:,k], uhist[:,k], dt)

    # Ak1, Bk1 = DT_lin_model(model,Xn[k],Un[k],dt)
    # xhist[:,k+1] = Xn[k+1] + Ak1*(xhist[:,k] - Xn[k]) + Bk1*(uhist[:,k] - Un[k])
    
    # update variables: shift and copy the last step
    Xprev[end] = Xhorizon[k][end]
    Xprev[1:end-1] = Xhorizon[k]
    Uprev[1:end-1] = Uhorizon[k][2:end]
    Uprev[end] = Uhorizon[k][end]
end

# Plot solution
Usln = hcat(Usln...)
Xsln = hcat(Xsln...)
plot_scene("OSQP-online-obs")

print("\nOSQP SIMULATION FINISHED!!!\n")
##
include("visualization.jl")
vis = rover_vis_init()
##
rover_vis_run(vis, Xsln, xhist, Xhorizon, xc, yc, rc, obs_flag, xf, 
                visextralines=true)
