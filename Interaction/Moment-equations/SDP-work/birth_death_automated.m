clc
% load data
data = load("data_BD.mat");
% extract capture scaling
B = data.B;
% extract orders
Nd = double(data.Nd);
Nda = double(data.Nda);
Ndbd = double(data.Ndbd);
% select sample
s = 2;
% extract bounds
key = sprintf("sample%d", s);
bounds = data.(key);
y_L = bounds(1, 1:Ndbd);
y_U = bounds(2, 1:Ndbd);
% moment vector
y = sdpvar(Nd, 1);
Cs = [y >= 0; y(1) == 1];
% moment bounds
Cs = [Cs, y_L(:) <= B * y(1:Ndbd)];
Cs = [Cs, B * y(1:Ndbd) <= y_U(:)];
% moment matrices
M0 = y(data.M0);
M1 = y(data.M1);
M2 = y(data.M2);
Cs = [Cs, M0 >= 0, M1 >= 0, M2 >= 0];
% rates
R = 5;
K = 100;
k = sdpvar(R, 1);
Cs = [Cs, k >= 0];
% Cs = [Cs, k(1) == 1];
Cs = [Cs, k <= K];
% moment equations
for alpha = 1:Nda
    key = sprintf("A%d", alpha - 1);
    A = data.(key);
    Cs = [Cs, k' * A * y == 0];
end
% optimize
sol = optimize(Cs, 1, sdpsettings('solver', 'bmibnb', 'bmibnb.maxtime', 10, 'bmibnb.maxiter', 1000))

if sol.problem == 0
    % Extract and display value
    solution = value(y)
    
else
    sol.info
    yalmiperror(sol.problem)
end
