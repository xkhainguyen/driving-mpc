"""OSQP with input constraint and obstacle constraint
Once the horizon hits the end of the trajectory, shrink the horizon. Try warm-start
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

print("\nSIMULATION START!!!")

# Get reference trajectory
dt = 1/20   # 10 - 20Hz
N = 155
tf = (N-1)*dt

steer_rate_lim=0.3
accel_lim=3

model = Rover(ref=:rear)
n,m = size(model);
Nx, Nu = size(model);

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
# xc = SA[6, 10, 12]
# yc = SA[5, 12, 5]
# rc = SA[2, 2, 2]    #try 1

xc = SA[6, 9, 12]
yc = SA[5, 12, 5]
rc = SA[2, 2, 2] 
rc_plan = SA[2.1, 2.1, 2.1]    #try 1

cir = CircleConstraint(n, xc, yc, rc_plan)
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
Un = copy(Usln)
Xn = copy(Xsln)  # ref state
K, P = riccati_recursion(model, Xn, Un)

thist = Array(range(0,dt*(N-1), step=dt))

# Initialize OSQP
Nh = 30
Nh′ = Nh

global x_ws = zeros(Nh*(Nx+Nu))

xhist = zeros(n,N)
xhist[:,1] .= xs[:,1]
uhist = zeros(m,N-1)

# u_bnd1 = [Inf;Inf]
u_bnd1 = 1.2*u_bnd
#

for k = 1:(N-1)  
    global Nh′
    println(k)
    if (k <= N-Nh)
        ubar0 = mpc_controller_obs2(k, xhist[:,k], Nh, model, P, Xn, Un, dt, x_bnd, u_bnd1)
    else
        Nh′ = Nh′ - 1
        ubar0 = mpc_controller_obs2(k, xhist[:,k], Nh′, model, P, Xn, Un, dt, x_bnd, u_bnd1; warm_start=false)
    end
    uhist[:,k] .= Un[:,k] + ubar0
    
    # xhist[:,k+1] .= DT_non_model(model, xhist[:,k], uhist[:,k], dt)

    Ak1, Bk1 = DT_lin_model(model,Xn[:,k],Un[:,k],dt)
    xhist[:,k+1] .= Xn[:,k+1] + Ak1*(xhist[:,k] - Xn[:,k]) + Bk1*ubar0
end

# Plot solution
plot_scene("OSQP-ws-obs-shrink")

print("SIMULATION FINISHED!!!")
##
indx = 5
plot(thist, 10*Xn[indx,:])
plot!(thist, 10*xhist[indx,:])
##
include("visualization.jl")
vis = rover_vis_init()
##
rover_vis_run(vis, Xsln, xhist, Nh, xc, yc, rc, xf, visextralines=true)
