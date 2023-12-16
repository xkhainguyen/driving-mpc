"""
Generate reference trajectory for second-order bicycle model
In TinyMPC project `bicycle_tvlqr`
"""

using TrajectoryOptimization
import RobotZoo.Rover
using StaticArrays, LinearAlgebra
using Plots
using RobotDynamics
using Altro
const RD = RobotDynamics
using SparseArrays

include("utils.jl")
include("my_mpc.jl")

printstyled("\nTRAJOPT SIMULATION START!!!\n", color=:cyan, bold=true)

dt = 1/10   # 10 - 20Hz
tf = 10.
N = Int(tf/dt) + 1 

steer_rate_lim = 0.5
accel_lim = 2

model = Rover(ref=:rear)
Nx, Nu = size(model); n,m = size(model);

x0 = @SVector zeros(Nx)
xf = @SVector [20, 20, π/4, 0, 0]

# Objective function
Q = 1e-1*Diagonal(@SVector ones(n)) 
Qf = 1e4*Diagonal(@SVector ones(n))
R = 1e1*Diagonal(@SVector ones(m)) 
obj = LQRObjective(Q,R,Qf,xf,N);

# Create Empty ConstraintList
conSet = ConstraintList(n,m,N)

# Circle constraints
xc = SA[14, 10, 16]
yc = SA[6, 10, 15]
rc = SA[1.5, 1.5, 1.5] 
rc_plan = rc*1.4   #try 1
cir = CircleConstraint(n, xc, yc, rc_plan)
add_constraint!(conSet, cir, 2:N-1)

# Control-state Bounds
u_bnd = [accel_lim; steer_rate_lim]
x_bnd = [Inf;Inf;Inf;4;0.7]
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
set_options!(solver, show_summary=true)
initial_controls!(solver, U0)
solve!(solver);

Usln = Vector.(controls(solver))
Xsln = Vector.(states(solver))
  
Un = deepcopy(Usln)
Xn = deepcopy(Xsln)  # ref state

xhist = zeros(n,N)
x0_off = [0;0.0;0;0;0]           # initial offset
xhist[:,1] = x0 + x0_off
xhist2 = deepcopy(xhist)
uhist = zeros(m,N-1)

Usln = hcat(Usln...)
Xsln = hcat(Xsln...)

# Plot solution
Plots.CURRENT_PLOT.nullableplot = nothing
layout = @layout [a; b]
p1 = plot(Xsln[1,1:end], Xsln[2,1:end], color="orange",linewidth=4,label="ref",xticks=x0[1]-1:xf[1]+1,yticks=x0[2]-1:xf[2]+1,legend=:bottomright,aspect_ratio=:equal, xlabel="x position", ylabel="y position", xlims=(x0[1]-1,xf[1]+1), ylims=(x0[2]-1,xf[2]+1))    

scatter!([x0[1]],[x0[2]],color = "blue", label = "start", markershape=:rect, markersize = 10,alpha = 0.8)
scatter!([xf[1]],[xf[2]],color = "orange", label = "goal", markershape=:star, 
markersize = 12,alpha = 0.8)
circles = circle.([xc[1:end], yc[1:end], rc[1:end]]...)
plot!(circles, c=:gray, label="") 

# Plot control inputs
time = dt*(0:N-2)
p2 = plot(time, Usln[1,1:end], label="accel",lw=8,linestyle=:dot,
color = "dark cyan")
plot!(time, Usln[2,1:end], label="steer_rate", xlabel="time (s)", 
lw=8,linestyle=:dot, color = "Chocolate")

# Plot rates
time = dt*(0:N-1)
p3 = plot(time, Xsln[4,1:end], label="vel",lw=8,linestyle=:dot,
color = "dark cyan")
plot!(time, Xsln[5,1:end], label="steer",lw=8, linestyle=:dot, 
color = "Chocolate")
# plot!(time, x_bnd[5], color = "orange")
display(plot(p1, p3, layout=layout, size=(500,1100)))

print("OSQP SIMULATION FINISHED!!!")

## Save reference trajectory to files
open("xref_data.txt","a") do io
    for line in Xsln
        println(io, line)
    end
end

open("uref_data.txt","a") do io
    for line in Usln
        println(io, line)
    end
end