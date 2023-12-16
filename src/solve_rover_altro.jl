"""Altro-based NMPC with input constraint and obstacle constraint
- Redesigned version
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

print("\nALTRO SIMULATION START!!!")

# Get reference trajectory
dt = 1/5   # 10 - 20Hz
N = 141
Nh = 40
tf = (N-1)*dt

steer_rate_lim=0.3
accel_lim=3

model = Rover(ref=:rear)
Nx, Nu = size(model);

x0 = @SVector zeros(Nx)
xf = @SVector [20, 20, π/4, 0, 0] 

# Objective function
Q = 1e-1*Diagonal(@SVector ones(Nx)) 
Qf = 1e2*Diagonal(@SVector ones(Nx))
R = 1e1*Diagonal(@SVector ones(Nu)) 
obj = LQRObjective(Q,R,Qf,xf,Nh);

# Create Empty ConstraintList
conSet = ConstraintList(Nx,Nu,Nh)

# Circle constraints
xc = SA[14, 10, 16]
yc = SA[6, 10, 15]
rc = SA[1.5, 1.5, 1.5] 
rc_plan = SA[2.1, 2.1, 2.1]    #try 1

cir = CircleConstraint(Nx, xc, yc, rc_plan)
add_constraint!(conSet, cir, 2:Nh-1)

# Control-state Bounds
u_bnd = [accel_lim; steer_rate_lim]
x_bnd = [Inf;Inf;Inf;40;π/3]
bnd = BoundConstraint(Nx, Nu, x_min=-x_bnd, x_max=x_bnd, 
                    u_min=-u_bnd, u_max=u_bnd)
add_constraint!(conSet, bnd, 1:Nh-1)
# Goal Constraint
# goal = GoalConstraint(xf, 1:2)
# add_constraint!(conSet, goal, Nh)

prob = Problem(model, obj, x0, tf, xf=xf, constraints=conSet)

u0 = @SVector fill(0.0, Nu)
U0 = [u0 for k = 1:Nh-1]
initial_controls!(prob, U0)
rollout!(prob);

solver = ALTROSolver(prob)
set_options!(solver, show_summary=false)
initial_controls!(solver, U0)
solve!(solver);
X_altro = states(solver)
U_altro = controls(solver)
Usln = hcat(Vector.(U_altro)...)
Xsln = hcat(Vector.(X_altro)...)

# Rover simulator
thist = Array(range(0,dt*(N-1), step=dt))
xhist = zeros(Nx, N)
uhist = zeros(Nu, N-1)
x0_off = [0;0;0;0;0]
xhist[:,1] .= x0 + x0_off

# Run MPC
for k = 1:N-1
    println(k)
    TrajectoryOptimization.set_initial_state!(prob, xhist[:,k])
    solver = ALTROSolver(prob)
    # prob.x0 .= xhist[:,k]
    initial_controls!(solver, U_altro)
    initial_states!(solver, X_altro)
    set_options!(solver, reset_duals=false, penalty_initial=1.0)
    solve!(solver)
    # @show X_altro .= states(solver)
    U_altro .= controls(solver)
    uhist[:,k] = U_altro[1]
    xhist[:,k+1] .= DT_non_model(model, xhist[:,k], uhist[:,k], dt)
end

# Plot solution
Plots.CURRENT_PLOT.nullableplot = nothing
layout = @layout [a; b]

p1 = plot!(xhist[1,1:end], xhist[2,1:end], color="green",linewidth=3,label=name,xticks=x0[1]-1:xf[1]+1,yticks=x0[2]-1:xf[2]+1,legend=:bottomright,aspect_ratio=:equal, xlabel="x position", ylabel="y position", xlims=(x0[1]-1,xf[1]+1), ylims=(x0[2]-1,xf[2]+1))    

scatter!([x0[1]],[x0[2]],color = "blue", label = "start", markershape=:rect, markersize = 10,alpha = 0.8)
scatter!([xf[1]],[xf[2]],color = "green", label = "goal", markershape=:star, 
markersize = 12,alpha = 0.8)

circles = circle.([xc[1:end], yc[1:end], rc[1:end]]...)
plot!(circles, c=:red, label="") 

# Plot control inputs
time = dt*(0:N-2)
p2 = plot(time, uhist[1,1:end], label="accel+",lw=4, color = "cyan")
plot!(time, 10*uhist[2,1:end], label="10*steer+", xlabel="time (s)", 
ylabel="control input", lw=4, reuse=false, color = "orange")

display(plot(p1, p2, layout=layout, size=(500,1100)))

print("OSQP SIMULATION FINISHED!!!")

##
include("visualization.jl")
vis = rover_vis_init()
##
rover_vis_run(vis, Xsln, xhist, Nh, xc, yc, rc, xf, visextralines=true)
