# using Revise
using TrajectoryOptimization
import RobotZoo.Rover
using StaticArrays, LinearAlgebra
using Plots
using RobotDynamics
using Altro
const RD = RobotDynamics
# noise (accel) to model (RK for loop) TV-LQR
# reuse solution

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

# Circle constraints
xc = SA[5, 10, 12]
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
set_options!(solver, show_summary=true)
initial_controls!(solver, U0)
solve!(solver);

Usln = hcat(Vector.(controls(solver))...)
Xsln = hcat(Vector.(states(solver))...)

# Rover simulator
xs = zeros(5,N)
xs[:,1] .= x0
# dmodel = RD.DiscretizedDynamics{RD.RK4}(model)
# for i in 1:N-1
#     # u = 3*rand(2) .- 0.5   
#     xs[:,i+1]  = RD.discrete_dynamics(dmodel, xs[:,i]+0.01*randn(5), Usln[:,i], 0, dt)
# end

# Plot solution
layout = @layout [a; b]

p1 = plot(xs[1,1:end], xs[2,1:end], legend=false, xlabel="x position", ylabel="y position", color="gray",markershape=:circle, markersize = 5,widen=1.1)
plot!(Xsln[1,1:end], Xsln[2,1:end], color="green")
scatter!([x0[1]],[x0[2]],color = "yellow", label = "", markershape=:rect, markersize = 10,alpha = 0.8)
scatter!([xf[1]],[xf[2]],color = "green", label = "", markershape=:diamond, markersize = 10,alpha = 0.8)
scatter!(xc, yc, markershape=:circle, color ="red",markersize=rc*15)

# Plot control inputs
time = gettimes(solver)[1:end-1]

p2 = plot(time, Usln[1,1:end], label="accel",lw=3)
plot!(time, 10*Usln[2,1:end], label="10*steer", xlabel="time (s)", ylabel="control input", lw=3, reuse=false)

plot(p1, p2, layout=layout, size=(600,600))
