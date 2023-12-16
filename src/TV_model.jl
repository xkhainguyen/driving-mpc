function A(x)
    A = [0 0 -x[5]*sin(x[3]) 0 cos(x[3]);
        0 0 x[5]*cos(x[3]) 0 sin(x[3]);
        0 0 0 x[5]/(cos(x[4])^2) tan(x[4]);
        zeros(2, 5)]
end
function B()
    B = [zeros(3,2);
        0 1;
        1 0]
end

# x̄ = x - xn
# ū = u - un
# dx̄ = A(xn)*x̄ + B()*ū

