p = sdpvar(2, 2)
Fr = [p(1)*p(3) == p(2); p >= 0] %[p(1)*p(2) >= 2; p >= 0; p <= 10]
optimize(Fr, 1, sdpsettings('solver', 'bmibnb'))
value(p)