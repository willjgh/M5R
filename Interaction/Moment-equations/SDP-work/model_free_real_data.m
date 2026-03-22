clc
% load data
data = load("real_bounds.mat");
% extract capture scaling
B = data.B;
% select sample
s = 8;
key = sprintf("sample%d", s);
bounds = data.(key);
% extract bounds
y_L = bounds(1, 1:6);
y_U = bounds(2, 1:6);
% moment vector
y = sdpvar(6, 1);
Cs = [y(1) == 1];
% moment bounds
Cs = [Cs, y_L(:) <= B * y(:)];
Cs = [Cs, B * y(:) <= y_U(:)];
% moment matrices
M0 = [y(1) y(2) y(3); y(2) y(4) y(5); y(3) y(5) y(6)];
M1 = [y(2)];
M2 = [y(3)];
Cs = [Cs, M0 >= 0, M1 >= 0, M2 >= 0];
% factorization
Cs = [Cs, y(5) == y(2) * y(3)];
% optimize
sol = optimize(Cs, 1, sdpsettings('solver', 'bmibnb'))

if sol.problem == 0
    % Extract and display value
    solution = value(y)
    
else
    sol.info
    yalmiperror(sol.problem)
end
