function DT_non_model(model,xk,uk,h)
    """
    Discrete-time nonlinear model
    4th order runge-kutta method
    """
    k1 = RD.dynamics(model, xk, uk)
    k2 = RD.dynamics(model, xk+0.5*h*k1, uk)
    k3 = RD.dynamics(model, xk+0.5*h*k2, uk)
    k4 = RD.dynamics(model, xk+h*k3, uk)
    x = xk + h*(k1+2k2+2k3+k4)/6
    return x
end


function DT_lin_model(model,xk,uk,h)
    """
    Discrete-time linear model
    """
    A = Altro.ForwardDiff.jacobian(x -> DT_non_model(model,x,uk,h), xk)
    B = Altro.ForwardDiff.jacobian(u -> DT_non_model(model,xk,u,h), uk)
    return A, B
end


function CT_lin_model(model, x, u)
    """
    Continuous-time linear model
    """
    # At = [0 0 -x[4]*sin(x[3]) cos(x[3])     0;
    # 0 0  x[4]*cos(x[3]) sin(x[3])   0;
    # 0 0           0 tan(x[5]) x[4]*(tan(x[5])^2 + 1);
    # 0 0           0      0                 0;
    # 0 0           0       0                  0]
    # Bt = [zeros(3,2);1 0; 0 1]
    At = Altro.ForwardDiff.jacobian(xd -> RD.dynamics(model,xd,u), x)
    Bt = Altro.ForwardDiff.jacobian(ud -> RD.dynamics(model,x,ud), u)

    return At, Bt
end


function simulator(model, xk, uk, ulim, h)
    """
    High-fidelity simulator
    """
    for i in eachindex(uk)
        if uk[i] > ulim[i]
            uk[i] = ulim[i]
        elseif uk[i] < -ulim[i]
            uk[i] = -ulim[i]
        end
    end

    k1 = RD.dynamics(model, xk, uk)
    k2 = RD.dynamics(model, xk+0.5*h*k1, uk)
    k3 = RD.dynamics(model, xk+0.5*h*k2, uk)
    k4 = RD.dynamics(model, xk+h*k3, uk)
    x = xk + h*(k1+2k2+2k3+k4)/6
    return x
end

# Create a circle where (x,y) is the center and r the radius
function circle(x, y, r=1; n=30)
    θ = 0:360÷n:360
    Plots.Shape(r*sind.(θ) .+ x, r*cosd.(θ) .+ y)
end

function plot_scene(name)
    # Plot solution
    Plots.CURRENT_PLOT.nullableplot = nothing
    layout = @layout [a; b]
    # p1 = plot(xs[1,1:end], xs[2,1:end], color="gray",linewidth=3,label="open-loop")
    p1 = plot!(Xsln[1,1:end], Xsln[2,1:end], color="blue",linewidth = 5,linestyle=:dot,label="nominal")

    plot!(xhist[1,1:end], xhist[2,1:end], color="orange",linewidth=4,label=name,xticks=x0[1]-1:xf[1]+1,yticks=x0[2]-1:xf[2]+1,legend=:bottomright,aspect_ratio=:equal, xlabel="x position", ylabel="y position", xlims=(x0[1]-1,xf[1]+1), ylims=(x0[2]-1,xf[2]+1))    
    
    scatter!([x0[1]],[x0[2]],color = "blue", label = "start", markershape=:rect, markersize = 10,alpha = 0.8)
    scatter!([xf[1]],[xf[2]],color = "orange", label = "goal", markershape=:star, 
    markersize = 12,alpha = 0.8)
    plot!(xhist2[1,1:end], xhist2[2,1:end]) # plot nonlinear
    circles = circle.([xc[1:end], yc[1:end], rc[1:end]]...)
    plot!(circles, c=:gray, label="") 
    
    # Plot control inputs
    time = dt*(0:N-2)
    p2 = plot(time, Usln[1,1:end], label="accel",lw=8,linestyle=:dot,
    color = "dark cyan")
    plot!(time, uhist[1,1:end], label="accel+",lw=4, color = "cyan")
    plot!(time, 10*Usln[2,1:end], label="10*steer", xlabel="time (s)", 
    lw=8,linestyle=:dot, color = "Chocolate")
    plot!(time, 10*uhist[2,1:end], label="10*steer+", xlabel="time (s)", 
    ylabel="control input", lw=4, reuse=false, color = "orange")

    # Plot rates
    time = dt*(0:N-1)
    p3 = plot(time, Xsln[4,1:end], label="vel",lw=8,linestyle=:dot,
    color = "dark cyan")
    plot!(time, xhist[4,1:end], label="vel+",lw=4, color = "cyan")
    plot!(time, Xsln[5,1:end], label="steer",lw=8, linestyle=:dot, 
    color = "Chocolate")
    plot!(time, xhist[5,1:end], label="steer+", xlabel="time (s)", 
    ylabel="control input", lw=4, reuse=false, color = "orange")
    # plot!(time, x_bnd[5], color = "orange")
    display(plot(p1, p3, layout=layout, size=(500,1100)))
end

function plot_scene2(name)
    # Plot solution
    Plots.CURRENT_PLOT.nullableplot = nothing
    layout = @layout [a; b]
    # p1 = plot(xs[1,1:end], xs[2,1:end], color="gray",linewidth=3,label="open-loop")
    p1 = plot!(Xsln[1,1:end], Xsln[2,1:end], color="blue",linewidth = 5,linestyle=:dot,label="nominal")

    plot!(xhist[1,1:end], xhist[2,1:end], color="green",linewidth=3,label=name,xticks=x0[1]-1:xf[1]+1,yticks=x0[2]-1:xf[2]+1,legend=:bottomright,aspect_ratio=:equal, xlabel="x position", ylabel="y position", xlims=(x0[1]-1,xf[1]+1), ylims=(x0[2]-1,xf[2]+1))    
    
    scatter!([x0[1]],[x0[2]],color = "blue", label = "start", markershape=:rect, markersize = 10,alpha = 0.8)
    scatter!([xf[1]],[xf[2]],color = "green", label = "goal", markershape=:star, 
    markersize = 12,alpha = 0.8)
    
    circles = circle.([xc[1:end], yc[1:end], rc[1:end]]...)
    plot!(circles, c=:red, label="") 
    
    # Plot control inputs
    time = dt*(0:N-2)
    p2 = plot(time, Usln[1,1:end], label="accel",lw=8,linestyle=:dot,
    color = "dark cyan")
    plot!(time, uhist[1,1:end], label="accel+",lw=4, color = "cyan")
    plot!(time, 10*Usln[2,1:end], label="10*steer", xlabel="time (s)", 
    lw=8,linestyle=:dot, color = "Chocolate")
    plot!(time, 10*uhist[2,1:end], label="10*steer+", xlabel="time (s)", 
    ylabel="control input", lw=4, reuse=false, color = "orange")

    display(plot(p1, p2, layout=layout, size=(500,1100)))
end