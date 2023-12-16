function riccati_recursion(model, Xn, Un)
    """Riccati recursion

    :param model:
    :param Xn:
    :param Un:
    :return K:
    :return P:
    """
    Nx, Nu = size(model)

    P = zeros(Nx, Nx, N)
    K = zeros(Nu, Nx, N - 1)
    P[:, :, N] .= Qf

    for k = (N-1):-1:1
        AA, BB = DT_lin_model(model, Xn[ k], Un[ k], dt)
        K[:, :, k] = (R + BB' * P[:, :, k+1] * BB) \ (BB' * P[:, :, k+1] * AA)
        P[:, :, k] = Q + AA' * P[:, :, k+1] * (AA - BB * K[:, :, k])
    end

    return K, P
end

function check_obstacle(x)
    """Calculate current obstacle constraints
    Ax >= b

    :param xk: current full state 
    :return A: diag of gradient of fi w.r.t x
    :return b: vector of RHS
    """
    delta_d = 7 # 7

    No = length(xc)
    obs_active = zeros(Int8, 0)
    X = x[1:2]

    for k = 1:No
        Xc = [xc[k], yc[k]]
        d = sqrt((X[1] - Xc[1])^2 + (X[2] - Xc[2])^2)
        # constraint is active if within `delta_d` from the obstacle boundary
        if d <= rc_plan[k] + delta_d
            append!(obs_active, floor(Int8, k))
            obs_flag[k] = true
        end
    end

    return obs_active
end

function add_obstacle(xref, obs_active)
    """Calculate future obstacle constraints
    AΔx >= b

    :param xk: current full state 
    :return A: diag of gradient of fi w.r.t x
    :return b: vector of RHS
    """
    A = zeros(0, 2)     # empty size 0x2
    b = zeros(0)        # empty size 0x1
    Xref = xref[1:2]

    for k in obs_active
        Xc = [xc[k]; yc[k]]
        vecXC = Xc - Xref
        d = sqrt(vecXC[1]^2 + vecXC[2]^2)
        vecXI = vecXC * (d - rc_plan[k]) / d
        Xi = Xref + vecXI
        # f = d^2 - r^2; ∇f*(Xref + Δx - Xi) > 0
        ∇f = [2 * (Xi[1] - Xc[1]) 2 * (Xi[2] - Xc[2])]
        A = vcat(A, ∇f)
        # @show Xi - Xref
        # @show ∇f*(Xi - Xref)
        b = vcat(b, ∇f * (Xi - Xref))
    end

    return A, b
end

function add_obstacle2(xref, obs_active)
    """Calculate future obstacle constraints
    AΔx >= b

    :param xk: current full state 
    :return A: diag of gradient of fi w.r.t x
    :return b: vector of RHS
    """
    A = zeros(0, 2)     # empty size 0x2
    b = zeros(0)        # empty size 0x1
    Xref = xref[1:2]

    for k in obs_active
        Xc = [xc[k]; yc[k]]
        vecXC = Xc - Xref
        d = sqrt(vecXC[1]^2 + vecXC[2]^2)
        vecXI = vecXC * (d - rc_plan[k]) / d
        Xi = Xref + vecXI
        # f = d^2 - r^2; ∇f*(Xref + Δx - Xi) > 0
        ∇f = [2 * (Xi[1] - Xc[1]) 2 * (Xi[2] - Xc[2])]
        A = vcat(A, ∇f)
        # @show Xi - Xref
        # @show ∇f*(Xi - Xref)
        b = vcat(b, ∇f * Xi)
    end

    return A, b
end

function mpc_controller_obs1(k, x, Nh, model, Xn, Un, dt, x_bnd, u_bnd;
    warm_start=true)
    """Run MPC controller on current state
    Constraints include obstacle and input

    :param k: current time step
    :param x: current state (initial for MPC)
    :param Nh: horizon length (time step)
    :param model: Rover model
    :param P: cost to go
    :param Xn: reference states
    :param Un: reference controls
    :param dt: time interval
    :param x_bnd: state bounds
    :param u_bnd: control bounds
    
    :return: xbar0: first control
    """

    Nx, Nu = size(model)
    No = 0
    knot = 1
    # Cost weights
    Q = sparse(1.0 * I(Nx))
    R = sparse(1.0 * I(Nu))

    # Build QP problem
    U_all = kron(I(Nh), [I zeros(Nu, Nx)]) # Matrix that picks out all u from decision variable
    H = kron(I(Nh), blockdiag(R, Q))
    b = zeros(Nh * (Nx + Nu))               # linear cost vector, dont need here

    if mod(k, knot) == 0
        @show obs_active = check_obstacle(x)
        No = length(obs_active)
    end
    # Constraints
    lb = zeros(Nh * (Nx + Nu + No))
    ub = zeros(Nh * (Nx + Nu + No))

    Ak1, Bk1 = DT_lin_model(model, Xn[ k], Un[ k], dt)
    lb[1:Nx] = -Ak1 * (x - Xn[ k])
    ub[1:Nx] = -Ak1 * (x - Xn[ k]) # equality constraints => dynamics
    C = kron(I(Nh), [Bk1 -I(n)])
    # control limits
    lb[size(C, 1).+(1:Nu)] = -u_bnd - Un[ k]
    ub[size(C, 1).+(1:Nu)] = u_bnd - Un[ k]
    # obstacles
    Xo = kron(I(Nh), zeros(No, Nx + Nu))
    if mod(k, knot) == 0
        A_obs, b_obs = add_obstacle(Xn[ k+1], obs_active)
        ub[size(C, 1)+size(U_all, 1).+(1:No)] = Inf * ones(No)
        lb[size(C, 1)+size(U_all, 1).+(1:No)] = b_obs
        Xo[1:No, 1:(Nx+Nu)] = A_obs * [0 0 1 0 0 0 0
            0 0 0 1 0 0 0]
    end

    for j = 1:Nh-1
        Ak, Bk = DT_lin_model(model, Xn[ k+j], Un[ k+j], dt)
        C[(j*Nx).+(1:Nx), (j*(Nx+Nu)-Nx).+(1:Nx+Nu)] = [Ak Bk]
        lb[j*Nu+size(C, 1).+(1:Nu)] = -u_bnd - Un[ k+j]
        ub[j*Nu+size(C, 1).+(1:Nu)] = u_bnd - Un[ k+j]
        # obstacles
        if mod(k, knot) == 0
            A_obs, b_obs = add_obstacle(Xn[ k+j+1], obs_active)
            ub[j*No+size(C, 1)+size(U_all, 1).+(1:No)] = Inf * ones(No)
            lb[j*No+size(C, 1)+size(U_all, 1).+(1:No)] = b_obs
            Xo[(j*No).+(1:No), (j*(Nx+Nu)).+(1:(Nx+Nu))] = A_obs * [0 0 1 0 0 0 0
                0 0 0 1 0 0 0]
        end
    end

    D = [C; U_all; Xo]   # constraints: [dynamics; control; x obstacles]
    # PH = sparse(P[:,:,k+Nh])
    # H[end-Nx+1:end, end-Nx+1:end] = PH

    prob = OSQP.Model()
    OSQP.setup!(prob; P=H, q=b, A=D, l=lb, u=ub, verbose=false, eps_abs=1e-4, eps_rel=1e-4, polish=1, warm_start=warm_start)

    if warm_start == true
        OSQP.warm_start_x!(prob, x_ws)
    end
    #Solve QP
    results = OSQP.solve!(prob)
    if results.info.status != :Solved
        @show results.info.status
    end

    xhorizon = kron(I(Nh), [zeros(Nx, Nu) I(Nx)]) * results.x
    xhorizon = reshape(xhorizon, (Nx, Nh)) + Xn[ k+1:k+Nh]
    xhorizon = [xhorizon[:, i] for i in 1:size(xhorizon, 2)]

    uhorizon = kron(I(Nh), [I(Nu) zeros(Nu, Nx)]) * results.x
    uhorizon = reshape(uhorizon, (Nu, Nh)) + Un[ k+1:k+Nh]
    uhorizon = [uhorizon[:, i] for i in 1:size(uhorizon, 2)]
    global x_ws = vcat(results.x[Nx+Nu+1:end], results.x[end-Nx-Nu+1:end])
    return xhorizon, uhorizon
end

function mpc_controller_obs2(k, x, Nh, model, Xn, Un, dt, x_bnd, u_bnd;
    warm_start=true)
    """Run MPC controller on current state
    Constraints include obstacle (state), rate and input
    """

    Nx, Nu = size(model)
    No = 0
    knot = 1
    # Cost weights
    Q = sparse(1.0 * I(Nx))
    R = sparse(1.0 * I(Nu))

    # Build QP problem
    U_all = kron(I(Nh), [I zeros(Nu, Nx)]) # Matrix that picks out all u from decision variable
    X_45 = kron(I(Nh), [zeros(2, Nu) [0 0 0 1 0; 0 0 0 0 1]])   # Picks out all x4 and x5 
    H = kron(I(Nh), blockdiag(R, Q))
    b = zeros(Nh * (Nx + Nu))   # linear cost vector, dont need here

    if mod(k, knot) == 0
        obs_active = check_obstacle(x)
        No = length(obs_active)
    end
    # All constraints
    lb = zeros(Nh * (Nx + Nu + 2 + No))
    ub = zeros(Nh * (Nx + Nu + 2 + No))

    Ak1, Bk1 = DT_lin_model(model, Xn[k], Un[k], dt)
    lb[1:Nx] = -Ak1 * (x - Xn[ k])
    ub[1:Nx] = -Ak1 * (x - Xn[ k]) # equality constraints => dynamics
    C = kron(I(Nh), [Bk1 -I(n)])
    # control limits
    lb[size(C, 1).+(1:Nu)] = -u_bnd - Un[ k]
    ub[size(C, 1).+(1:Nu)] = u_bnd - Un[ k]
    # rate limit
    lb[size(C, 1)+size(U_all, 1).+(1:2)] = -x_bnd[4:5] - Xn[k][4:5]
    ub[size(C, 1)+size(U_all, 1).+(1:2)] = x_bnd[4:5] - Xn[k][4:5]  
    # obstacles
    Xo = kron(I(Nh), zeros(No, Nx + Nu))
    if mod(k, knot) == 0
        A_obs, b_obs = add_obstacle(Xn[ k+1], obs_active)
        ub[size(C, 1)+size(X_45, 1)+size(U_all, 1).+(1:No)] = Inf * ones(No)
        lb[size(C, 1)+size(X_45, 1)+size(U_all, 1).+(1:No)] = b_obs
        Xo[1:No, 1:(Nx+Nu)] = A_obs * [0 0 1 0 0 0 0
            0 0 0 1 0 0 0]
    end

    for j = 1:Nh-1
        Ak, Bk = DT_lin_model(model, Xn[ k+j], Un[ k+j], dt)
        C[(j*Nx).+(1:Nx), (j*(Nx+Nu)-Nx).+(1:Nx+Nu)] = [Ak Bk]
        lb[j*Nu+size(C, 1).+(1:Nu)] = -u_bnd - Un[k+j]
        ub[j*Nu+size(C, 1).+(1:Nu)] = u_bnd - Un[k+j]
        # rate limit  
        lb[j*2+size(C, 1)+size(U_all, 1).+(1:2)] = -x_bnd[4:5] - Xn[k+j][4:5]
        ub[j*2+size(C, 1)+size(U_all, 1).+(1:2)] = x_bnd[4:5] - Xn[k+j][4:5]
        # obstacles
        if mod(k, knot) == 0
            A_obs, b_obs = add_obstacle(Xn[ k+j+1], obs_active)
            ub[j*No+size(C, 1)+size(X_45, 1)+size(U_all, 1).+(1:No)] = Inf * ones(No)
            lb[j*No+size(C, 1)+size(X_45, 1)+size(U_all, 1).+(1:No)] = b_obs
            Xo[(j*No).+(1:No), (j*(Nx+Nu)).+(1:(Nx+Nu))] = A_obs * [0 0 1 0 0 0 0
                0 0 0 1 0 0 0]
        end
    end

    D = [C; U_all; X_45; Xo]   # constraints: [dynamics; control; x obstacles]
   
    prob = OSQP.Model()
    OSQP.setup!(prob; P=H, q=b, A=D, l=lb, u=ub, verbose=false, eps_abs=1e-4, eps_rel=1e-4, polish=1, warm_start=warm_start)

    if warm_start == true
        OSQP.warm_start_x!(prob, x_ws)
    end
    #Solve QP
    results = OSQP.solve!(prob)
    if results.info.status != :Solved
        @show results.info.status
    end

    xhorizon = kron(I(Nh), [zeros(Nx, Nu) I(Nx)]) * results.x
    xhorizon = reshape(xhorizon, (Nx, Nh)) 
    xhorizon = [xhorizon[:, i] + Xn[k+i] for i in 1:size(xhorizon, 2)]

    uhorizon = kron(I(Nh), [I(Nu) zeros(Nu, Nx)]) * results.x
    uhorizon = reshape(uhorizon, (Nu, Nh)) 
    uhorizon = [uhorizon[:, i] + Un[k+i] for i in 1:size(uhorizon, 2)]
    global x_ws = vcat(results.x[Nx+Nu+1:end], results.x[end-Nx-Nu+1:end])
    return xhorizon, uhorizon
end

function mpc_controller_obs3(k, x, Nh, model, Xn, Un, dt, x_bnd, u_bnd;
    warm_start=true)
    """Run MPC controller on current state
    Constraints include obstacle (state), rate and input
    """
    #TODO: CONVERT OBS HARD TO SOFT CONSTRAINTS, 
    #TODO: ENFORCE ANGLE CONSTRAINT?
    Nx, Nu = size(model)
    No = 0
    knot = 1
    # Cost weights
    Q = sparse(1.0 * I(Nx))
    R = sparse(1.0 * I(Nu))

    # Build QP problem
    U_all = kron(I(Nh), [I zeros(Nu, Nx)]) # Matrix that picks out all u from decision variable
    X_345 = kron(I(Nh), [zeros(3, Nu) [0 0 1 0 0; 
                                        0 0 0 1 0;
                                        0 0 0 0 1]]) 
    H = kron(I(Nh), blockdiag(R, Q))
    b = zeros(Nh * (Nx + Nu))   # linear cost vector, dont need here

    if mod(k, knot) == 0
        obs_active = check_obstacle(x)
        No = length(obs_active)
    end
    # All constraints
    lb = zeros(Nh * (Nx + Nu + 3 + No))
    ub = zeros(Nh * (Nx + Nu + 3 + No))

    Ak1, Bk1 = DT_lin_model(model, Xn[ k], Un[ k], dt)
    lb[1:Nx] = -Ak1 * (x - Xn[ k])
    ub[1:Nx] = -Ak1 * (x - Xn[ k]) # equality constraints => dynamics
    C = kron(I(Nh), [Bk1 -I(n)])
    # control limits
    lb[size(C, 1).+(1:Nu)] = -u_bnd - Un[ k]
    ub[size(C, 1).+(1:Nu)] = u_bnd - Un[ k]
    # rate limit
    theta_bnd = 0.174*3
    theta_ubnd = theta_bnd
    theta_lbnd = -theta_bnd
    lb[size(C, 1)+size(U_all, 1).+(1:3)] = vcat([theta_lbnd; -x_bnd[4:5] - Xn[4:5,k]])
    ub[size(C, 1)+size(U_all, 1).+(1:3)] = vcat([theta_ubnd; x_bnd[4:5] - Xn[4:5,k]])    
    # obstacles
    Xo = kron(I(Nh), zeros(No, Nx + Nu))
    if mod(k, knot) == 0
        A_obs, b_obs = add_obstacle(Xn[ k+1], obs_active)
        ub[size(C, 1)+size(X_345, 1)+size(U_all, 1).+(1:No)] = Inf * ones(No)
        lb[size(C, 1)+size(X_345, 1)+size(U_all, 1).+(1:No)] = b_obs
        Xo[1:No, 1:(Nx+Nu)] = A_obs * [0 0 1 0 0 0 0
            0 0 0 1 0 0 0]
    end

    for j = 1:Nh-1
        Ak, Bk = DT_lin_model(model, Xn[ k+j], Un[ k+j], dt)
        C[(j*Nx).+(1:Nx), (j*(Nx+Nu)-Nx).+(1:Nx+Nu)] = [Ak Bk]
        lb[j*Nu+size(C, 1).+(1:Nu)] = -u_bnd - Un[ k+j]
        ub[j*Nu+size(C, 1).+(1:Nu)] = u_bnd - Un[ k+j]
        # rate limit  
        theta_ubnd = theta_bnd
        theta_lbnd = -theta_bnd
        lb[j*3+size(C, 1)+size(U_all, 1).+(1:3)] = vcat([theta_lbnd;-x_bnd[4:5] - Xn[4:5,k+j]])
        ub[j*3+size(C, 1)+size(U_all, 1).+(1:3)] = vcat([theta_ubnd;x_bnd[4:5] - Xn[4:5,k+j]])
        # obstacles
        if mod(k, knot) == 0
            A_obs, b_obs = add_obstacle(Xn[ k+j+1], obs_active)
            ub[j*No+size(C, 1)+size(X_345, 1)+size(U_all, 1).+(1:No)] = Inf * ones(No)
            lb[j*No+size(C, 1)+size(X_345, 1)+size(U_all, 1).+(1:No)] = b_obs
            Xo[(j*No).+(1:No), (j*(Nx+Nu)).+(1:(Nx+Nu))] = A_obs * [0 0 1 0 0 0 0
                0 0 0 1 0 0 0]
        end
    end

    D = [C; U_all; X_345; Xo]   # constraints: [dynamics; control; x obstacles]
   
    prob = OSQP.Model()
    OSQP.setup!(prob; P=H, q=b, A=D, l=lb, u=ub, verbose=false, eps_abs=1e-4, eps_rel=1e-4, polish=1, warm_start=warm_start)

    if warm_start == true
        OSQP.warm_start_x!(prob, x_ws)
    end
    #Solve QP
    results = OSQP.solve!(prob)
    if results.info.status != :Solved
        @show results.info.status
    end

    xhorizon = kron(I(Nh), [zeros(Nx, Nu) I(Nx)]) * results.x
    xhorizon = reshape(xhorizon, (Nx, Nh)) 
    xhorizon = [xhorizon[:, i] + Xn[k+i] for i in 1:size(xhorizon, 2)]

    uhorizon = kron(I(Nh), [I(Nu) zeros(Nu, Nx)]) * results.x
    uhorizon = reshape(uhorizon, (Nu, Nh)) 
    uhorizon = [uhorizon[:, i] + Un[k+i] for i in 1:size(uhorizon, 2)]
    global x_ws = vcat(results.x[Nx+Nu+1:end], results.x[end-Nx-Nu+1:end])
    return xhorizon, uhorizon
end

function mpc_controller_abs(k, x, Nh, model, Xn, Un, dt, x_bnd, u_bnd; warm_start=true)
    """Decision variables are X and U (absolute)
    No obstacle constraints
    """
    Nx, Nu = size(model)
    No = 0      # number of active obstacles
    knot = 1    # collocation point to enforce obs constraints
    
    # Cost weights
    Q = sparse(1e0 * I(Nx))
    R = sparse(1e0 * I(Nu))

    # Build QP problem
    H = kron(I(Nh), blockdiag(R, Q))     # Hessian of the cost
    q = kron(ones(Nh, 1), [-R * Un[k]; -Q * Xn[k+1]])    # linear cost
    for j=1:Nh-1
        q[j*(Nx+Nu).+(1:Nu+Nx)] = [-R * Un[k+j]; -Q * Xn[k+1+j]]
    end
    q = vec(q)

    # Matrix that picks out all u from decision variable
    U_all = kron(I(Nh), [I zeros(Nu, Nx)])
    # Picks out all x4 and x5 
    X_45 = kron(I(Nh), [zeros(2, Nu) [0 0 0 1 0; 0 0 0 0 1]])  

    if mod(k, knot) == 0
        obs_active = check_obstacle(x)
        No = length(obs_active)
    end

    # Constraint 
    lb = zeros(Nh * (Nx + Nu + 2 + No))
    ub = zeros(Nh * (Nx + Nu + 2 + No))
    # linearize about current state
    Ak1, Bk1 = DT_lin_model(model, Xn[k], Un[k], dt)
    # display(norm(x - Xn[k]))
    b_temp = -Ak1 * (x - Xn[k]) - Xn[k+1] + Bk1 * Un[k]
    lb[1:Nx] = b_temp
    ub[1:Nx] = b_temp #  equality constraints => dynamics
    C = kron(I(Nh), [Bk1 -I(n)])
    # control limits
    lb[size(C, 1).+(1:Nu)] = -u_bnd
    ub[size(C, 1).+(1:Nu)] = u_bnd
    # rate limit
    lb[size(C, 1)+size(U_all, 1).+(1:2)] = -x_bnd[4:5] 
    ub[size(C, 1)+size(U_all, 1).+(1:2)] = x_bnd[4:5] 
    # obstacles
    Xo = kron(I(Nh), zeros(No, Nx + Nu))

    if mod(k, knot) == 0
        A_obs, b_obs = add_obstacle2(Xn[k+1], obs_active)
        ub[size(C, 1)+size(X_45, 1)+size(U_all, 1).+(1:No)] = Inf * ones(No)
        lb[size(C, 1)+size(X_45, 1)+size(U_all, 1).+(1:No)] = b_obs
        Xo[1:No, 1:(Nx+Nu)] = A_obs * [0 0 1 0 0 0 0; 0 0 0 1 0 0 0]
    end

    for j = 1:Nh-1
        Ak, Bk = DT_lin_model(model, Xn[k+1+j], Un[k+j], dt)
        C[(j*Nx).+(1:Nx), (j*(Nx+Nu)-Nx).+(1:Nx+Nu)] = [Ak Bk]
        b_temp = -Xn[k+2+j] + Ak * Xn[k+1+j] + Bk * Un[k+j]
        lb[j*Nx.+(1:Nx)] = b_temp
        ub[j*Nx.+(1:Nx)] = b_temp #  equality constraints => dynamics 
        # control limit   
        lb[j*Nu+size(C, 1).+(1:Nu)] = -u_bnd
        ub[j*Nu+size(C, 1).+(1:Nu)] = u_bnd
        # rate limit  
        lb[j*2+size(C, 1)+size(U_all, 1).+(1:2)] = -x_bnd[4:5] 
        ub[j*2+size(C, 1)+size(U_all, 1).+(1:2)] = x_bnd[4:5] 
        # obstacles
        if mod(k, knot) == 0
            A_obs, b_obs = add_obstacle2(Xn[k+2+j], obs_active)
            ub[j*No+size(C, 1)+size(X_45, 1)+size(U_all, 1).+(1:No)] = Inf * ones(No)
            lb[j*No+size(C, 1)+size(X_45, 1)+size(U_all, 1).+(1:No)] = b_obs
            Xo[(j*No).+(1:No), (j*(Nx+Nu)).+(1:(Nx+Nu))] = A_obs * [0 0 1 0 0 0 0; 0 0 0 1 0 0 0]
        end
    end

    # constraints: [dynamics; input; rates; obstacles]
    D = [C; U_all; X_45; Xo]

    prob = OSQP.Model()
    OSQP.setup!(prob; P=H, q=q, A=D, l=lb, u=ub, verbose=false, eps_abs=1e-4, eps_rel=1e-4, polish=1, warm_start=warm_start)

    if warm_start == true
        # OSQP.warm_start!(prob; x=x_ws)
        OSQP.warm_start_x!(prob, x_ws)
    end

    #Solve QP
    results = OSQP.solve!(prob)
    if results.info.status != :Solved
        @show results.info.status
    end

    xhorizon, uhorizon = extract_states_controls(results, Nh, Nx, Nu)
    global x_ws = vcat(results.x[Nx+Nu+1:end], results.x[end-Nx-Nu+1:end])
    return xhorizon, uhorizon
end

function extract_states_controls(results, Nh, Nx, Nu)
    """Extract states and controls from OSQP results
    return `vector{vector}` type
    """
    # Extract states and convert to vector of vectors
    x = kron(I(Nh), [zeros(Nx, Nu) I(Nx)]) * results.x
    x = reshape(x, (Nx, Nh))
    x = [x[:, i] for i in 1:size(x, 2)]
    # Extract control and convert to vector of vectors
    u = kron(I(Nh), [I(Nu) zeros(Nu, Nx)]) * results.x
    u = reshape(u, (Nu, Nh))
    u = [u[:, i] for i in 1:size(u, 2)]
    return x, u
end

function mpc_controller_online(k, x, xf, Nh, model, Xprev, Uprev, dt, x_bnd, u_bnd; warm_start=true)
    """Decision variables are X and U (absolute)
    No obstacle constraints
    """
    Nx, Nu = size(model)
    No = 0      # number of active obstacles
    knot = Inf    # collocation point to enforce obs constraints

    # Cost weights
    Q = sparse(1e0 * I(Nx))
    # Q = sparse(Diagonal(1e0*[1, 1, 0, 0, 0]))
    R = sparse(1e0 * I(Nu))
    # Qf = Diagonal(100e-1*[1, 1, 0, 0, 0])
    Qf = sparse(2 * Q)
    uf = [0.0; 0]
    # Build QP problem
    H = kron(I(Nh), blockdiag(R, Q))     # Hessian of the cost
    H[end-Nx+1:end, end-Nx+1:end] = Qf     # terminal cost

    # Matrix that picks out all x from decision variable
    q = kron(ones(Nh, 1), [-R * uf; -Q * xf])    # linear cost
    q[end-Nx+1:end] = -Qf * xf                # terminal linear cost
    q = vec(q)

    # Matrix that picks out all u from decision variable
    U_all = kron(I(Nh), [I zeros(Nu, Nx)])
    # Picks out all x4 and x5 
    X_45 = kron(I(Nh), [zeros(2, Nu) [0 0 0 1 0; 0 0 0 0 1]])  

    if mod(k, knot) == 0
        @show obs_active = check_obstacle(x)
        No = length(obs_active)
    end

    # Constraint 
    lb = zeros(Nh * (Nx + Nu + 2 + No))
    ub = zeros(Nh * (Nx + Nu + 2 + No))
    # linearize about current state
    Ak1, Bk1 = DT_lin_model(model, Xprev[1], Uprev[1], dt)
    display(norm(x - Xprev[1]))
    b_temp = -Ak1 * (x - Xprev[1]) - Xprev[2] + Bk1 * Uprev[1]
    lb[1:Nx] .= b_temp
    ub[1:Nx] .= b_temp #  equality constraints => dynamics
    C = kron(I(Nh), [Bk1 -I(n)])
    # control limits
    lb[size(C, 1).+(1:Nu)] .= -u_bnd
    ub[size(C, 1).+(1:Nu)] .= u_bnd
    # rate limit
    lb[size(C, 1)+size(U_all, 1).+(1:2)] .= -x_bnd[4:5] 
    ub[size(C, 1)+size(U_all, 1).+(1:2)] .= x_bnd[4:5] 
    # obstacles
    Xo = kron(I(Nh), zeros(No, Nx + Nu))

    if mod(k, knot) == 0
        A_obs, b_obs = add_obstacle2(Xprev[2], obs_active)
        ub[size(C, 1)+size(X_45, 1)+size(U_all, 1).+(1:No)] = Inf * ones(No)
        lb[size(C, 1)+size(X_45, 1)+size(U_all, 1).+(1:No)] = b_obs
        Xo[1:No, 1:(Nx+Nu)] = A_obs * [0 0 1 0 0 0 0; 0 0 0 1 0 0 0]
    end

    for j = 1:Nh-1
        Ak, Bk = DT_lin_model(model, Xprev[1+j], Uprev[1+j], dt)
        C[(j*Nx).+(1:Nx), (j*(Nx+Nu)-Nx).+(1:Nx+Nu)] = [Ak Bk]
        b_temp = -Xprev[j+2] + Ak * Xprev[j+1] + Bk * Uprev[j+1]
        lb[j*Nx.+(1:Nx)] .= b_temp
        ub[j*Nx.+(1:Nx)] .= b_temp #  equality constraints => dynamics    
        lb[j*Nu+size(C, 1).+(1:Nu)] .= -u_bnd
        ub[j*Nu+size(C, 1).+(1:Nu)] .= u_bnd
        # rate limit  
        lb[j*2+size(C, 1)+size(U_all, 1).+(1:2)] .= -x_bnd[4:5] 
        ub[j*2+size(C, 1)+size(U_all, 1).+(1:2)] .= x_bnd[4:5] 
        # obstacles
        if mod(k, knot) == 0
            A_obs, b_obs = add_obstacle2(Xprev[j+2], obs_active)
            ub[j*No+size(C, 1)+size(X_45, 1)+size(U_all, 1).+(1:No)] = Inf * ones(No)
            lb[j*No+size(C, 1)+size(X_45, 1)+size(U_all, 1).+(1:No)] .= b_obs
            Xo[(j*No).+(1:No), (j*(Nx+Nu)).+(1:(Nx+Nu))] = A_obs * [0 0 1 0 0 0 0; 0 0 0 1 0 0 0]
        end
    end

    # constraints: [dynamics; input; rates; obstacles]
    D = [C; U_all; X_45; Xo]

    prob = OSQP.Model()
    OSQP.setup!(prob; P=H, q=q, A=D, l=lb, u=ub, verbose=false, eps_abs=1e-4, eps_rel=1e-4, polish=1, warm_start=warm_start)

    if warm_start == true
        # OSQP.warm_start!(prob; x=x_ws)
        # OSQP.warm_start_x!(prob, x_ws)
    end

    #Solve QP
    results = OSQP.solve!(prob)
    if results.info.status != :Solved
        @show results.info.status
    end

    xhorizon, uhorizon = extract_states_controls(results, Nh, Nx, Nu)
    global x_ws = vcat(results.x[Nx+Nu+1:end], results.x[end-Nx-Nu+1:end])
    return xhorizon, uhorizon
end

function mpc_controller_online2(k, x, xf, Nh, model, Xprev, Uprev, dt, x_bnd, u_bnd; warm_start=true)
    """Decision variables are X and U (absolute)
    No obstacle constraints
    Constrains small angle approximation
    """
    Nx, Nu = size(model)
    No = 0      # number of active obstacles
    knot = Inf    # collocation point to enforce obs constraints

    # Cost weights
    Q = sparse(1e0 * I(Nx))
    # Q = sparse(Diagonal(1e0*[1, 1, 0, 0, 0]))
    R = sparse(1e0 * I(Nu))
    # Qf = Diagonal(100e-1*[1, 1, 0, 0, 0])
    Qf = sparse(10 * Q)
    uf = [0.0; 0]
    # Build QP problem
    H = kron(I(Nh), blockdiag(R, Q))     # Hessian of the cost
    H[end-Nx+1:end, end-Nx+1:end] = Qf     # terminal cost

    # Matrix that picks out all x from decision variable
    q = kron(ones(Nh, 1), [-R * uf; -Q * xf])    # linear cost
    q[end-Nx+1:end] = -Qf * xf                # terminal linear cost
    q = vec(q)

    # Matrix that picks out all u from decision variable
    U_all = kron(I(Nh), [I zeros(Nu, Nx)])
    # Picks out all x4 and x5 
    X_345 = kron(I(Nh), [zeros(3, Nu) [0 0 1 0 0; 
                                        0 0 0 1 0;
                                        0 0 0 0 1]])  

    if mod(k, knot) == 0
        @show obs_active = check_obstacle(x)
        No = length(obs_active)
    end

    # Constraint 
    lb = zeros(Nh * (Nx + Nu + 3 + No))
    ub = zeros(Nh * (Nx + Nu + 3 + No))
    # linearize about current state
    Ak1, Bk1 = DT_lin_model(model, Xprev[1], Uprev[1], dt)
    display(norm(x - Xprev[1]))
    b_temp = -Ak1 * (x - Xprev[1]) - Xprev[2] + Bk1 * Uprev[1]
    lb[1:Nx] .= b_temp
    ub[1:Nx] .= b_temp #  equality constraints => dynamics
    C = kron(I(Nh), [Bk1 -I(n)])
    # control limits
    lb[size(C, 1).+(1:Nu)] .= -u_bnd
    ub[size(C, 1).+(1:Nu)] .= u_bnd
    # rate limit
    theta_ubnd = Xprev[2][3] + 0.174/2
    theta_lbnd = Xprev[2][3] - 0.174/2
    lb[size(C, 1)+size(U_all, 1).+(1:3)] .= vcat([theta_lbnd; -x_bnd[4:5]])
    ub[size(C, 1)+size(U_all, 1).+(1:3)] .= vcat([theta_ubnd; x_bnd[4:5]])
    # obstacles
    Xo = kron(I(Nh), zeros(No, Nx + Nu))

    if mod(k, knot) == 0
        A_obs, b_obs = add_obstacle2(Xprev[2], obs_active)
        ub[size(C, 1)+size(X_345, 1)+size(U_all, 1).+(1:No)] = Inf * ones(No)
        lb[size(C, 1)+size(X_345, 1)+size(U_all, 1).+(1:No)] = b_obs
        Xo[1:No, 1:(Nx+Nu)] = A_obs * [0 0 1 0 0 0 0; 0 0 0 1 0 0 0]
    end

    for j = 1:Nh-1
        Ak, Bk = DT_lin_model(model, Xprev[1+j], Uprev[1+j], dt)
        C[(j*Nx).+(1:Nx), (j*(Nx+Nu)-Nx).+(1:Nx+Nu)] = [Ak Bk]
        b_temp = -Xprev[j+2] + Ak * Xprev[j+1] + Bk * Uprev[j+1]
        lb[j*Nx.+(1:Nx)] .= b_temp
        ub[j*Nx.+(1:Nx)] .= b_temp #  equality constraints => dynamics    
        lb[j*Nu+size(C, 1).+(1:Nu)] .= -u_bnd
        ub[j*Nu+size(C, 1).+(1:Nu)] .= u_bnd
        # rate limit  
        theta_ubnd = Xprev[2+j][3] + 0.174/2
        theta_lbnd = Xprev[2+j][3] - 0.174/2
        lb[j*3+size(C, 1)+size(U_all, 1).+(1:3)] .= vcat([theta_lbnd; -x_bnd[4:5]])
        ub[j*3+size(C, 1)+size(U_all, 1).+(1:3)] .= vcat([theta_ubnd; x_bnd[4:5]])
        # obstacles
        if mod(k, knot) == 0
            A_obs, b_obs = add_obstacle2(Xprev[j+2], obs_active)
            ub[j*No+size(C, 1)+size(X_345, 1)+size(U_all, 1).+(1:No)] = Inf * ones(No)
            lb[j*No+size(C, 1)+size(X_345, 1)+size(U_all, 1).+(1:No)] .= b_obs
            Xo[(j*No).+(1:No), (j*(Nx+Nu)).+(1:(Nx+Nu))] = A_obs * [0 0 1 0 0 0 0; 0 0 0 1 0 0 0]
        end
    end

    # constraints: [dynamics; input; rates; obstacles]
    D = [C; U_all; X_345; Xo]

    prob = OSQP.Model()
    OSQP.setup!(prob; P=H, q=q, A=D, l=lb, u=ub, verbose=false, eps_abs=1e-4, eps_rel=1e-4, polish=1, warm_start=warm_start)

    if warm_start == true
        # OSQP.warm_start!(prob; x=x_ws)
        OSQP.warm_start_x!(prob, x_ws)
    end

    #Solve QP
    results = OSQP.solve!(prob)
    if results.info.status != :Solved
        @show results.info.status
    end

    xhorizon, uhorizon = extract_states_controls(results, Nh, Nx, Nu)
    global x_ws = vcat(results.x[Nx+Nu+1:end], results.x[end-Nx-Nu+1:end])
    return xhorizon, uhorizon
end