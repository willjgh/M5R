clc
% moment vector
y = sdpvar(6, 1)
Cs = [y(1) == 1]
% moment bounds
y_L = [1 2 2 6 4 6] % correctly becomes infeasible if E[X1X2] != E[X1]E[X2] in bounds
y_U = [1 2 2 6 4 6]
Cs = [Cs, y_L(:) <= y(:)]
Cs = [Cs, y(:) <= y_U(:)]
% moment matrices
M0 = [y(1) y(2) y(3); y(2) y(4) y(5); y(3) y(5) y(6)]
M1 = [y(2)]
M2 = [y(3)]
Cs = [Cs, M0 >= 0, M1 >= 0, M2 >= 0]
% factorization
Cs = [Cs, y(5) == y(2) * y(3)]
% optimize
sol = optimize(Cs, 1, sdpsettings('solver', 'bmibnb'))

if sol.problem == 0
    % Extract and display value
    solution = value(y)
    
else
    sol.info
    yalmiperror(sol.problem)
end
