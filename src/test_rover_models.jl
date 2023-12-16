# using Revise
using TrajectoryOptimization
import RobotZoo.Rover
using StaticArrays, LinearAlgebra
using Plots
using RobotDynamics
using Altro
const RD = RobotDynamics
using SparseArrays

include("visualization.jl")
include("utils.jl")

## SOLVER
dt = 1/20   # 10 - 20Hz
N = 151
tf = (N-1)*dt

steer_rate_lim=0.3
accel_lim=3

model = Rover(ref=:rear)
n,m = size(model);

x0 = @SVector zeros(n)
xf = @SVector [15, 15, π/4, 0, 0];  

# Objective function
Q = 1e-1*Diagonal(@SVector ones(n)) 
Qf = 1e3*Diagonal(@SVector ones(n))
R = 1e1*Diagonal(@SVector ones(m)) 
obj = LQRObjective(Q,R,Qf,xf,N);

# Create Empty ConstraintList
conSet = ConstraintList(n,m,N)

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
set_options!(solver, show_summary=true)
initial_controls!(solver, U0)
solve!(solver);

Usln = hcat(Vector.(controls(solver))...)
Xsln = hcat(Vector.(states(solver))...)

## Rover simulator
Xnonl = zeros(5,N)
Xnonl[:,1] .= Xsln[:,1]
Xlin = zeros(5,N)
Xlin[:,1] .= Xsln[:,1] - Xsln[:,1]
u = Usln*0.9
# dmodel2 = RD.DiscretizedDynamics{RD.RK4}(model2)
for k in 1:N-1
    Ak, Bk = DT_lin_model(model,Xsln[:,k],Usln[:,k],dt)
    Xlin[:,k+1] .= Xsln[:,k+1] + Ak*(Xlin[:,k]-Xsln[:,k]) + Bk*(u[:,k]-Usln[:,k])
    Xnonl[:,k+1] .= DT_non_model(model, Xnonl[:,k], u[:,k], dt)
end

# Plot solution

p1 = plot(Xnonl[1,1:end], Xnonl[2,1:end], xlabel="x position", ylabel="y position", color="gray",linewidth=3,widen=1.1,label="Nonlinear")
plot!(Xlin[1,1:end], Xlin[2,1:end], color="green",linewidth = 5,linestyle=:dot,label="Linear")
# plot!(Xsln[1,1:end], Xsln[2,1:end], color="blue",linewidth = 4,label="Linear")
display(p1)

print("SIMULATION FINISHED!!!")